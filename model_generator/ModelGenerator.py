#!/usr/bin/python
# -*- coding: utf-8 -*-

from SBMLParser import SBMLParser
from ModelCodeGenerator import ModelCodeGenerator
from constants import *
import os


class ModelStructure:

    def __init__(self, cid, parser, internal_reaction_sets=None, external_reaction_sets=None):
        """
        Defines the basic structure of a compartment. A compartment has at least one internal
        reaction set called the bulk (i.e. the compartment internal reactions).

        :param cid: The compartment id
        :param internal_reaction_sets: A list of all the compartment internal sets
        :param parser: A SBMLParser with the sbml file loaded
        :param external_reaction_sets:  A list of all the compartment external sets. Optional
        """

        if internal_reaction_sets is None:
            internal_reaction_sets = []

        if external_reaction_sets is None:
            external_reaction_sets = []

        self.id = cid
        self.space = {}
        self.routing_table = {}
        self.reaction_sets = {}

        # Compartment input-membrane mapping
        self.membrane_eic = {rsn: port
                             for rsn, port
                             in zip(internal_reaction_sets, range(internal_reaction_sets))}

        # Building routing table, each reaction can be reached using the port
        # that the routing table indicates. Each port communicates the space with
        # a different reaction_set
        #
        # Initial ports are for IC (internal reaction sets)
        internal_reaction_sets = [BULK] + internal_reaction_sets
        routing_sets = internal_reaction_sets
        port_numbers = range(len(routing_sets))
        self.routing_table = {(self.id, rsn): port
                              for rsn, port
                              in zip(routing_sets, port_numbers)}

        # Final ports are for EOC (external reaction sets from other compartments)
        port_numbers = range(len(self.routing_table),
                             len(self.routing_table) + len(external_reaction_sets))
        self.routing_table.update({external_cid: port
                                   for external_cid, port
                                   in zip(external_reaction_sets, port_numbers)})

        # Building reaction sets
        self.reaction_sets = {}
        for reaction_set in internal_reaction_sets:
            reaction_ids = parser.get_reaction_set_ids(cid, reaction_set)
            self.reaction_sets[reaction_set] = {rid: parser.reactions[rid] for rid in reaction_ids}

        # Building space
        reaction_locations = parser.get_reaction_locations(cid)
        self.space = {
            'species': parser.parse_compartments_species()[cid].keys(),
            'reaction_locations': reaction_locations,
            'enzymes': parser.get_enzymes(set(reaction_locations.keys()))
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

        special_comp_ids = [extra_cellular_id, periplasm_id, cytoplasm_id]
        cytoplasm_external_reaction_sets = [(periplasm_id, reaction_set)
                                            for reaction_set in [INNER, TRANS]]
        organelle_membranes = [(cid, MEMBRANE)
                               for cid in self.parser.get_compartments()
                               if cid not in special_comp_ids]
        cytoplasm_external_reaction_sets += organelle_membranes

        extra_cellular_external_reaction_sets = [(periplasm_id, reaction_set)
                                                 for reaction_set in [OUTER, TRANS]]

        self.periplasm = ModelStructure(periplasm_id, self.parser, [OUTER, INNER, TRANS])

        self.extra_cellular = ModelStructure(extra_cellular_id,
                                             self.parser,
                                             external_reaction_sets=extra_cellular_external_reaction_sets)

        self.cytoplasm = ModelStructure(cytoplasm_id,
                                        self.parser,
                                        external_reaction_sets=cytoplasm_external_reaction_sets)

        self.organelles = [ModelStructure(comp_id, self.parser, [MEMBRANE])
                           for comp_id in self.parser.get_compartments()
                           if comp_id not in special_comp_ids]

    def generate_reaction_sets(self, compartment):
        cid = compartment.id
        return {(cid, rsn): self.generate_reaction_set(cid, rsn, reaction_set)
                for rsn, reaction_set in compartment.reaction_sets.iteritems()}

    def generate_reaction_set(self, cid, rsn, reaction_set):
        reaction_set_id = cid + '_' + rsn

        router_routing_table = {rid: port
                                for rid, port
                                in zip(reaction_set.keys(), range(reaction_set))}
        router = self.coder.write_atomic_model(ROUTER_MODEL_CLASS,
                                               reaction_set_id,
                                               [self.parameter_file.name, reaction_set_id],
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

        return self.coder.write_coupled_model(reaction_set_id,
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
        sub_models = [space] + [reaction_set for (reaction_set, _, _) in reaction_sets.values()]

        # All membranes use output port 0 to send messages to the space, this is why the
        # coupled output port zero is not used, membranes uses port zero to send product to the
        # space.
        out_port_amount = max([out_port_amount
                               for (_, _, out_port_amount)
                               in reaction_sets.values()])
        in_ports = [(port, REACTANT_MESSAGE_TYPE, 'in')
                    for port in compartment.membrane_eic.values()]

        out_ports = [(port_number, PRODUCT_MESSAGE_TYPE, 'out')
                     for (_, _, output_port_amount) in reaction_sets.values()
                     for port_number in range(1, output_port_amount)]
        ports = in_ports + out_ports

        ic = []
        for (rs_cid, rs_name), port_number in compartment.routing_table.iteritems():
            # Organelle spaces do not send metabolites to other compartments without passing through
            # their membranes, that is the assert
            assert rs_cid == compartment.id
            ic += [(space, port_number, reaction_sets[rs_name][0], 0),
                   (reaction_sets[rs_name][0], 0, space, 0)]

        # If the output port amount is equal to 1, then the model only sends messages to the space
        # and there is no external communication, thus, it must not be linked in the EOC.
        # Because the port 0 of each reaction set goes to the space, the output ports range is
        # [1, .. , output port amount)
        eoc = [(reaction_set, port_number, port_number)
               for (reaction_set, _, output_port_amount) in reaction_sets.values()
               for port_number in range(1, output_port_amount)]

        eic = [(port_number, reaction_set, 0)
               for reaction_set, port_number
               in compartment.membrane_eic.iteritems()]

        return self.coder.write_coupled_model(compartment.id,
                                              sub_models,
                                              ports,
                                              eic,
                                              eoc,
                                              ic)

    # bulk compartments are cytoplasm and extracellular space, they represent the space where the
    # cell and the organelles live
    def generate_bulk_compartment(self, compartment):
        cid = compartment.id
        reaction_sets = self.generate_reaction_sets(compartment)
        assert len(reaction_sets) == 1  # cytoplasm has no membranes

        bulk = reaction_sets[BULK][0]

        output_port_amount = max(compartment.routing_table.values()) + 1
        space = self.coder.write_atomic_model(SPACE_MODEL_CLASS,
                                              cid,
                                              [self.parameter_file.name, cid],
                                              output_port_amount,
                                              1,
                                              REACTANT_MESSAGE_TYPE,
                                              PRODUCT_MESSAGE_TYPE)

        sub_models = [space, bulk]

        in_ports = [(0, PRODUCT_MESSAGE_TYPE, 'in')]
        out_ports = [(port_number, REACTANT_MESSAGE_TYPE, 'out')
                     for (rs_cid, _), port_number in compartment.routing_table.iteritems()
                     if rs_cid != cid]
        ports = in_ports + out_ports

        bulk_port_number = compartment.routing_table[(cid, BULK)]
        ic = [(space, bulk_port_number, bulk, 0), (bulk, 0, space, 0)]

        eoc = [(space, port_number, port_number)
               for (rs_cid, _), port_number in compartment.routing_table.iteritems()
               if rs_cid != cid]

        eic = [(0, space, 0)]

        return self.coder.write_coupled_model(cid,
                                              sub_models,
                                              ports,
                                              eic,
                                              eoc,
                                              ic)

    def generate_top(self):

        cytoplasm_model = self.generate_bulk_compartment(self.cytoplasm)
        extra_cellular_model = self.generate_bulk_compartment(self.extra_cellular)
        periplasm_model = self.generate_organelle_compartment(self.periplasm)

        sub_models = [cytoplasm_model, extra_cellular_model, periplasm_model]

        organelle_models = {}
        for organelle in self.organelles:
            organelle_models[organelle] = self.generate_organelle_compartment(organelle)
            sub_models.append(organelle_models[organelle])

        ports = []
        eoc = []
        eic = []
        ic = []
        for organelle in self.organelles:
            (organelle_model, _, output_port_amount) = organelle_models[organelle]

            cytoplasm_to_organelle_port = next(port_number
                                               for (cid, _), port_number
                                               in self.cytoplasm.routing_table.iteritems()
                                               if cid == organelle)
            ic.append([(organelle_model, 0, cytoplasm_model, 0),
                       (cytoplasm_model, cytoplasm_to_organelle_port, organelle_model, 0)])

            if output_port_amount > 1:
                extra_cellular_to_organelle = next(port_number
                                                   for (cid, _), port_number
                                                   in self.extra_cellular.routing_table.iteritems()
                                                   if cid == organelle)
                ic.append([(organelle_model, 1, extra_cellular_model, 0),
                           (extra_cellular_model, extra_cellular_to_organelle, organelle_model, 0)])

        return self.coder.write_coupled_model('cell', sub_models, ports, eic, eoc, ic)
