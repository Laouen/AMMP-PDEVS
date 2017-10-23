#!/usr/bin/python
# -*- coding: utf-8 -*-

from SBMLParser import SBMLParser
from ModelCodeGenerator import ModelCodeGenerator
from constants import *
import os


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
                port_number += 1
                self.routing_table[(external_cid, reaction_set)] = port_number

        # Building reaction sets
        self.reaction_sets = {}
        for reaction_set in internal_reaction_sets:
            reaction_ids = parser.get_reaction_set_ids(cid, reaction_set)
            self.reaction_sets[reaction_set] = {rid: parser.reactions[rid] for rid in reaction_ids}

        # Building space
        reaction_parameters = parser.get_reaction_parameters(cid)
        metabolites = {specie: parser.metabolite_amounts[specie]
                       for specie in parser.parse_compartments_species()[cid].keys()}
        self.space = {
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
                 parameters_path=os.curdir,
                 reactions=None,
                 enzymes=None):
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
        :param parameters_path: The path to the xml file where to write the model parameter values
        :param reactions: The parser.reactions loaded objects. Optional, used to avoid re parsing
        :param enzymes: T parser.enzymes luaded objects. Optional, used to avoid re parsing
        """

        self.parameter_file = open(parameters_path + os.sep + 'parameters.xml', 'w')

        self.coder = ModelCodeGenerator()
        self.parser = SBMLParser(sbml_file,
                                 extra_cellular_id,
                                 periplasm_id,
                                 cytoplasm_id,
                                 reactions,
                                 enzymes)

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
        rs_id = '_'.join([cid, rsn])

        router_routing_table = {rid: port
                                for rid, port
                                in zip(reaction_set.keys(), range(reaction_set))}
        router = self.coder.write_atomic_model(ROUTER_MODEL_CLASS,
                                               rs_id,
                                               [self.parameter_file.name, rs_id],
                                               len(router_routing_table),
                                               1,
                                               REACTANT_MESSAGE_TYPE,
                                               REACTANT_MESSAGE_TYPE)
        sub_models = [router]
        eic = [(0, router, 0)]
        eoc = []
        ic = []
        total_out_ports = 0
        for rid, parameters in reaction_set.iteritems():
            output_port_amount = max(parameters['routing_table'].values()) + 1
            total_out_ports = max(output_port_amount, total_out_ports)
            reaction = self.coder.write_atomic_model(REACTION_MODEL_CLASS,
                                                     rid,
                                                     [self.parameter_file.name, rid],
                                                     output_port_amount,
                                                     1,
                                                     PRODUCT_MESSAGE_TYPE,
                                                     REACTANT_MESSAGE_TYPE)
            ic.append((router, router_routing_table[rid], reaction, 0))
            eoc += [(reaction, port_number, port_number)
                    for port_number in range(output_port_amount)]
            sub_models.append(reaction)

        ports = [(0, REACTANT_MESSAGE_TYPE, 'in')]
        ports += [(port_number, PRODUCT_MESSAGE_TYPE, 'out')
                  for port_number in range(total_out_ports)]

        return self.coder.write_coupled_model(rs_id,
                                              sub_models,
                                              ports,
                                              eic,
                                              eoc,
                                              ic)

    def generate_organelle_compartment(self, compartment):

        reaction_sets = self.generate_reaction_sets(compartment)

        space = self.coder.write_atomic_model(SPACE_MODEL_CLASS,
                                              compartment.id,
                                              [self.parameter_file.name, compartment.id],
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
            ic += [(space, port_number, rs_model_name, 0), (rs_model_name, 0, space, 0)]

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
                eoc.append((rs_model_name, port_number, port_number))
                if port_number not in out_port_numbers:
                    out_port_numbers.append(port_number)
                    out_ports.append((port_number, PRODUCT_MESSAGE_TYPE, 'out'))

        in_ports = []
        eic = []
        for rs_name, port_number in compartment.membrane_eic.iteritems():
            rs_model_name = reaction_sets[(compartment.id, rs_name)][0]
            eic.append((port_number, rs_model_name, 0))
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

        bulk = reaction_sets[BULK][0]

        # port numbers starts in zero -> port amount = max{port numbers} + 1.
        output_port_amount = max(compartment.routing_table.values()) + 1
        space = self.coder.write_atomic_model(SPACE_MODEL_CLASS,
                                              cid,
                                              [self.parameter_file.name, cid],
                                              output_port_amount,
                                              1,
                                              REACTANT_MESSAGE_TYPE,
                                              PRODUCT_MESSAGE_TYPE)

        sub_models = [space, bulk]

        bulk_port_number = compartment.routing_table[(cid, BULK)]
        ic = [(space, bulk_port_number, bulk, 0), (bulk, 0, space, 0)]

        out_ports = []
        eoc = []
        for (rs_cid, _), port_number in compartment.routing_table.iteritems():
            if rs_cid != cid:
                eoc.append((space, port_number, port_number))
                out_ports.append((port_number, REACTANT_MESSAGE_TYPE, 'out'))

        in_ports = [(0, PRODUCT_MESSAGE_TYPE, 'in')]
        eic = [(0, space, 0)]

        return self.coder.write_coupled_model(cid,
                                              sub_models,
                                              in_ports + out_ports,
                                              eic,
                                              eoc,
                                              ic)

    def generate_top(self):

        cytoplasm_model = self.generate_bulk_compartment(self.cytoplasm)[0]
        extra_cellular_model = self.generate_bulk_compartment(self.extra_cellular)[0]
        periplasm_model = self.generate_organelle_compartment(self.periplasm)[0]

        sub_models = [cytoplasm_model, extra_cellular_model, periplasm_model]

        organelle_models = {}
        for cid, model_structure in self.organelles.iteritems():
            organelle_models[cid] = self.generate_organelle_compartment(model_structure)
            sub_models.append(organelle_models[cid][0])  # append model name

        ic = []
        for cid, (model_name, _, output_port_amount) in organelle_models.iteritems():
            for rsn, rs_port_number in self.organelles[cid].membrane_eic.iteritems():
                c_port_number = self.cytoplasm.routing_table[(cid, rsn)]
                ic.append((cytoplasm_model, c_port_number, model_name, rs_port_number))
                ic.append((model_name, 0, cytoplasm_model, 0))

                if output_port_amount > 1:
                    e_port_number = self.cytoplasm.routing_table[(cid, rsn)]
                    ic.append((extra_cellular_model, e_port_number, model_name, rs_port_number))
                    ic.append((model_name, 1, extra_cellular_model, 0))

        return self.coder.write_coupled_model('cell', sub_models, [], [], [], ic)
