#!/usr/bin/python
# -*- coding: utf-8 -*-

from SBMLParser import SBMLParser
from ModelCodeGenerator import ModelCodeGenerator
from DynamicModelCodeGenerator import DynamicModelCodeGenerator
from XMLParametersGenerator import XMLParametersGenerator
from constants import *
from itertools import islice


def chunks(data, SIZE=150):
    it = iter(data)
    for i in xrange(0, len(data), SIZE):
        yield {k: data[k] for k in islice(it, SIZE)}


class ModelStructure:

    def __init__(self, cid, parser, membranes=None, external_reaction_sets=None):
        """
        Defines the basic structure of a compartment. A compartment has at least one internal
        reaction set called the bulk (i.e. the compartment internal reactions).

        :param cid: The compartment id.
        :param membranes: A list of all the compartment internal reaction sets less the bulk.
        :param parser: A SBMLParser with the sbml file loaded.
        :param external_reaction_sets:  A dictionary of all the related external reaction sets
        grouped by their compartments. Optional.
        :type external_reaction_sets: dict[str, list[str]]
        """

        if membranes is None:
            membranes = []

        if external_reaction_sets is None:
            external_reaction_sets = {}

        self.id = cid
        self.space = {}
        self.routing_table = {}
        self.reaction_sets = {}

        # Compartment input-membrane mapping
        port_numbers = range(len(membranes))
        self.membrane_eic = {rsn: port
                             for rsn, port
                             in zip(membranes, port_numbers)}

        # Building routing table, each reaction can be reached using the port
        # that the routing table indicates. Each port communicates the space with
        # a different reaction_set
        #
        # The bulk reaction set is always present but is not a membrane, that is why it is not
        # considered in the input-membrane mapping.
        #
        # Initial ports are for IC (internal reaction sets)
        internal_reaction_sets = [BULK] + membranes
        port_numbers = range(len(internal_reaction_sets))
        self.routing_table = {(self.id, rsn): port
                              for rsn, port
                              in zip(internal_reaction_sets, port_numbers)}

        # Final ports are for EOC (external reaction sets from other compartments)
        port_number = len(self.routing_table)
        for external_cid, reaction_sets in external_reaction_sets.iteritems():
            for reaction_set in reaction_sets:
                self.routing_table[(external_cid, reaction_set)] = port_number
                port_number += 1

        # Building reaction sets
        self.reaction_sets = {}
        for reaction_set in internal_reaction_sets:
            reaction_ids = parser.get_reaction_set_rids(cid, reaction_set)
            self.reaction_sets[reaction_set] = {rid: parser.reactions[rid] for rid in reaction_ids}

        # Building space
        reaction_parameters = parser.get_reaction_parameters(cid)
        metabolites = {specie: parser.metabolite_amounts[specie]
                       for specie in parser.parse_compartments_species()[cid].keys()}
        self.space_parameters = {
            'cid': cid,
            'interval_time': parser.interval_times[cid],
            'metabolites': metabolites,
            'reaction_parameters': reaction_parameters,
            'enzymes': parser.get_enzymes(reaction_parameters.keys())
        }


class ModelGenerator:
    """
    Generates a model structure from a SBML file using the SBMLParser to get the information
    """

    def __init__(self,
                 sbml_file,
                 extra_cellular_id,
                 periplasm_id,
                 cytoplasm_id,
                 json_model=None,
                 groups_size=150):
        """
        Generates the whole model structure and generates a .cpp file with a cadmium model from the
        generated structure

        :param sbml_file: Path to the sbml file to parse
        :param extra_cellular_id: The extra cellular ID, this ID must be one of the compartments
        IDs from the sbml file
        :param periplasm_id: The periplasm ID, this ID must be one of the compartments IDs from
        the sbml file
        :param cytoplasm_id: The cytoplasm ID, this ID must be one of the compartments IDs from
        the sbml file
        :param json_model: The parser exported as json. Optional, used to avoid re parsing
        :param groups_size: The size of the reaction set groups
        """

        self.groups_size = groups_size
        self.parameter_writer = XMLParametersGenerator()
        self.coder = DynamicModelCodeGenerator()
        self.parser = SBMLParser(sbml_file,
                                 extra_cellular_id,
                                 periplasm_id,
                                 cytoplasm_id,
                                 json_model=json_model)

        special_compartment_ids = [extra_cellular_id, periplasm_id, cytoplasm_id]

        c_external_reaction_sets = {cid: [MEMBRANE]
                                    for cid in self.parser.get_compartments()
                                    if cid not in special_compartment_ids}
        c_external_reaction_sets[periplasm_id] = [INNER, TRANS]

        e_external_reaction_sets = {periplasm_id: [OUTER, TRANS]}

        # Compartment model structures
        self.periplasm = ModelStructure(periplasm_id,
                                        self.parser,
                                        membranes=[OUTER, INNER, TRANS])

        self.extra_cellular = ModelStructure(extra_cellular_id,
                                             self.parser,
                                             external_reaction_sets=e_external_reaction_sets)

        self.cytoplasm = ModelStructure(cytoplasm_id,
                                        self.parser,
                                        external_reaction_sets=c_external_reaction_sets)

        self.organelles = {comp_id: ModelStructure(comp_id, self.parser, [MEMBRANE])
                           for comp_id in self.parser.get_compartments()
                           if comp_id not in special_compartment_ids}

    def generate_reaction_sets(self, compartment):
        cid = compartment.id
        return {(cid, rsn): self.generate_reaction_set(cid, rsn, reaction_set)
                for rsn, reaction_set in compartment.reaction_sets.iteritems()}

    def generate_reaction_set(self, cid, rsn, reaction_set):
        reaction_set_id = '_'.join([cid, rsn])

        # Reaction sets are separated in groups, this is because the C++ problem compiling large tuples
        # Currently the maximum group amount is 150 and each group can hold 150 reactions, thus, each compartment
        # can have a maximum of 150 * 150 = 22500 reactions.
        groups = chunks(reaction_set, SIZE=self.groups_size)
        groups_reaction_ids = []
        routing_table = {}
        
        # Routers xml parameters
        port = 0
        for group in groups:
            reaction_ids = group.keys()
            groups_reaction_ids.append(reaction_ids)

            group_id = '_'.join([cid, rsn, str(port)])
            routing_table.update({rid: port for rid in reaction_ids})
            port += 1

            port_numbers = range(len(reaction_ids))
            group_routing_table = {rid: port_number for rid, port_number in zip(reaction_ids, port_numbers)}
            self.parameter_writer.add_router(group_id, group_routing_table)

        self.parameter_writer.add_router(reaction_set_id, routing_table)

        # Reaction xml parameters 
        for rid, parameters in reaction_set.iteritems():
            self.parameter_writer.add_reaction(rid, parameters)

        return self.coder.write_reaction_set(cid, rsn, groups_reaction_ids, self.parameter_writer.xml_file_path)

    def generate_organelle_compartment(self, compartment):

        reaction_sets = self.generate_reaction_sets(compartment)

        self.parameter_writer.add_space(compartment.id,
                                        compartment.space_parameters,
                                        compartment.routing_table)
        space = self.coder.write_atomic_model(SPACE_MODEL_CLASS,
                                              compartment.id,
                                              [self.parameter_writer.xml_file_path, compartment.id],
                                              len(reaction_sets),
                                              1,
                                              REACTANT_MESSAGE_TYPE,
                                              PRODUCT_MESSAGE_TYPE)

        sub_models = [space] + [rs_model_name for (rs_model_name, _, _) in reaction_sets.values()]

        ic = []
        for rs_address, port_number in compartment.routing_table.iteritems():
            # Organelle spaces do not send metabolites to other compartments without passing through
            # their membranes, that is the assert.
            assert rs_address[0] == compartment.id
            rs_model_name = reaction_sets[rs_address][0]
            ic += [(space, space, port_number, rs_model_name, 'pmgbp::models::reaction', 0), (rs_model_name, 'pmgbp::models::reaction', 0, space, space, 0)]

        # If the output port amount is equal to 1, then the model only sends messages to the space
        # and there is no communication with the extra cellular space, thus, it must not be linked
        # in the EOC.
        # Because the port 0 of each reaction set goes to the space, the output ports range is
        # [1, ... , output_port_amount)
        out_port_numbers = []
        out_ports = []
        eoc = []
        for (rs_model_name, _, output_port_amount) in reaction_sets.values():
            for port_number in range(1, output_port_amount):
                eoc.append((rs_model_name, 'pmgbp::models::reaction', port_number, port_number))
                if port_number not in out_port_numbers:
                    out_port_numbers.append(port_number)
                    out_ports.append((port_number, PRODUCT_MESSAGE_TYPE, 'out'))

        in_ports = []
        eic = []
        for rs_name, port_number in compartment.membrane_eic.iteritems():
            rs_model_name = reaction_sets[(compartment.id, rs_name)][0]
            eic.append((rs_model_name, 'pmgbp::models::reaction', 0, port_number))
            in_ports.append((port_number, REACTANT_MESSAGE_TYPE, 'in'))

        return self.coder.write_coupled_model(compartment.id,
                                              sub_models,
                                              out_ports + in_ports,
                                              eic,
                                              eoc,
                                              ic)

    # Bulk compartments are cytoplasm and extracellular space, they represent the space where the
    # cell and the organelles live
    def generate_bulk_compartment(self, compartment):
        cid = compartment.id
        reaction_sets = self.generate_reaction_sets(compartment)
        assert len(reaction_sets) == 1  # cytoplasm has no membranes

        bulk = reaction_sets[(cid, BULK)][0]

        self.parameter_writer.add_space(cid,
                                        compartment.space_parameters,
                                        compartment.routing_table)
        # port numbers starts in zero -> port amount = max{port numbers} + 1.
        output_port_amount = max(compartment.routing_table.values()) + 1
        space = self.coder.write_atomic_model(SPACE_MODEL_CLASS,
                                              cid,
                                              [self.parameter_writer.xml_file_path, cid],
                                              output_port_amount,
                                              1,
                                              REACTANT_MESSAGE_TYPE,
                                              PRODUCT_MESSAGE_TYPE)

        sub_models = [space, bulk]

        bulk_port_number = compartment.routing_table[(cid, BULK)]
        ic = [(space, space, bulk_port_number, bulk, 'pmgbp::models::reaction', 0), (bulk, 'pmgbp::models::reaction', 0, space, space, 0)]

        out_ports = []
        eoc = []
        for (rs_cid, _), port_number in compartment.routing_table.iteritems():
            if rs_cid != cid:
                eoc.append((space, space, port_number, port_number))
                out_ports.append((port_number, REACTANT_MESSAGE_TYPE, 'out'))

        in_ports = [(0, PRODUCT_MESSAGE_TYPE, 'in')]
        eic = [(space, space, 0, 0)]

        return self.coder.write_coupled_model(cid,
                                              sub_models,
                                              in_ports + out_ports,
                                              eic,
                                              eoc,
                                              ic)

    def generate_top(self, top='cell'):

        cytoplasm_model = self.generate_bulk_compartment(self.cytoplasm)[0]
        extra_cellular_model = self.generate_bulk_compartment(self.extra_cellular)[0]
        periplasm_model = self.generate_organelle_compartment(self.periplasm)[0]

        sub_models = [cytoplasm_model, extra_cellular_model, periplasm_model]

        organelle_models = {}
        for cid, model_structure in self.organelles.iteritems():
            organelle_models[cid] = self.generate_organelle_compartment(model_structure)
            sub_models.append(organelle_models[cid][0])  # append model name

        ic = self.generate_top_bulk_ic(cytoplasm_model,
                                       self.cytoplasm.id,
                                       self.cytoplasm.routing_table,
                                       periplasm_model)

        ic += self.generate_top_bulk_ic(extra_cellular_model,
                                        self.extra_cellular.id,
                                        self.extra_cellular.routing_table,
                                        periplasm_model)

        for cid, (model_name, _, output_port_amount) in organelle_models.iteritems():
            for rsn, rs_port_number in self.organelles[cid].membrane_eic.iteritems():
                c_port_number = self.cytoplasm.routing_table[(cid, rsn)]
                ic.append((cytoplasm_model, c_port_number, model_name, rs_port_number))
                # *1) organelle output port 1 always goes to cytoplasm
                ic.append((model_name, model_name, 1, cytoplasm_model, cytoplasm_model, 0))

                if output_port_amount > 1:
                    e_port_number = self.cytoplasm.routing_table[(cid, rsn)]
                    ic.append((extra_cellular_model, extra_cellular_model, e_port_number, model_name, model_name, rs_port_number))
                    # *1) organelle output port 2 always goes to extracellular
                    ic.append((model_name, model_name, 2, extra_cellular_model, extra_cellular_model, 0))

        self.parameter_writer.save_xml()
        cell_coupled = self.coder.write_coupled_model(top, sub_models, [], [], [], ic)

        return cell_coupled

    def generate_top_bulk_ic(self, bulk_model, bulk_cid, bulk_routing_table, periplasm_model):
        ic = []

        # this is the same logic as *1)
        periplasm_oport_number = 1 if bulk_cid == self.cytoplasm.id else 2
        ic.append((periplasm_model, periplasm_model, periplasm_oport_number, bulk_model, bulk_model, 0))

        for (cid, rsn), c_port_number in bulk_routing_table.iteritems():
            if cid == bulk_cid:
                continue

            if cid == self.periplasm.id:
                periplasm_port_number = self.periplasm.membrane_eic[rsn]
                ic.append((bulk_model, bulk_model, c_port_number, periplasm_model, periplasm_model, periplasm_port_number))
        return ic

    def end_model(self):
        self.coder.end_model()

    # WARNING: deprecated
    def generate_reaction_group(self, group_id, reaction_group):  
        # write parameters to the xml file
        reaction_ids = reaction_group.keys()
        port_numbers = range(len(reaction_ids))
        routing_table = {rid: port_number for rid, port_number in zip(reaction_ids, port_numbers)}
        self.parameter_writer.add_router('router_' + group_id, routing_table)

        for rid, parameters in reaction_group.iteritems():
            self.parameter_writer.add_reaction('reaction_' + rid, parameters)

        # NOTE: the reaction_id order is important because it matches the port number, thus, the reaction_ids order
        # must not be modified once the routing_table is declared.
        return self.coder.write_reaction_group_template(group_id, reaction_ids, self.parameter_writer.xml_file_path) 
