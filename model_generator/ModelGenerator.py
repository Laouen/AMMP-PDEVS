#!/usr/bin/python
# -*- coding: utf-8 -*-

from SBMLParser import SBMLParser
from ModelCodeGenerator import ModelCodeGenerator


class ModelStructure:

    def __init__(self, comp_id, parser, internal_reaction_sets=None, external_reaction_sets=None):
        """
        Defines the basic structure of a compartment. A compartment has at least one internal
        reaction set, the bulk (i.e. the compartment internal reactions)

        :param comp_id: The compartment id
        :param internal_reaction_sets: A list of all the compartment internal sets
        :param parser: A SBMLParser with the sbml file loaded
        :param external_reaction_sets:  A list of all the compartment external sets. Optional
        """

        if internal_reaction_sets is None:
            internal_reaction_sets = []

        if external_reaction_sets is None:
            external_reaction_sets = []

        self.id = comp_id
        self.space = {}
        self.routing_table = {}
        self.reaction_sets = {}

        # Building routing table, each reaction can be reached using the port
        # that the routing table indicates. Each port communicates the space with
        # a different reaction_set
        internal_reaction_sets.append('bulk')
        routing_sets = internal_reaction_sets + external_reaction_sets
        rs_len = range(len(routing_sets))
        self.routing_table = {(comp_id, r): p for r, p in zip(routing_sets, rs_len)}

        # Building reaction sets
        self.reaction_sets = {}
        for reaction_set in internal_reaction_sets:
            reaction_ids = parser.get_reaction_set_ids(comp_id, reaction_set)
            self.reaction_sets[reaction_set] = {rid: parser.reactions[rid] for rid in reaction_ids}

        # Building space
        reaction_locations = parser.get_reaction_locations(comp_id)
        self.space = {
            'species': parser.parse_compartments_species()[comp_id].keys(),
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
        :param reactions: The parser.reactions loaded objects. Optional, used to avoid re parsing
        :param enzymes: T parser.enzymes luaded objects. Optional, used to avoid re parsing
        """

        self.coder = ModelCodeGenerator()
        self.parser = SBMLParser(sbml_file,
                                 extra_cellular_id,
                                 periplasm_id,
                                 cytoplasm_id,
                                 reactions,
                                 enzymes)

        periplasm_sets = ['outer', 'inner', 'trans']
        special_comp_ids = [extra_cellular_id, periplasm_id, cytoplasm_id]

        # TODO: correctly set the cytoplasm and extra cellular space external reaction sets
        self.periplasm = ModelStructure(periplasm_id, self.parser, periplasm_sets)
        self.extra_cellular = ModelStructure(extra_cellular_id, self.parser)
        self.cytoplasm = ModelStructure(cytoplasm_id, self.parser)
        self.organelles = [ModelStructure(comp_id, self.parser, ['membrane'])
                           for comp_id in self.parser.get_compartments()
                           if comp_id not in special_comp_ids]

    def generate_reaction_sets(self, compartment):

        return [self.generate_reaction_set(compartment.id, rsn, reaction_set)
                for rsn, reaction_set in compartment.reaction_sets.iteritems()]

    def generate_reaction_set(self, cid, rsn, reaction_set):
        reaction_set_id = cid + '_' + rsn
        router = self.coder.write_atomic_model('router',
                                               reaction_set_id,
                                               [],  # TODO: put correct parameters
                                               len(reaction_set),
                                               1,
                                               'Reactant',
                                               'Reactant')
        submodels = [router]
        eic = [(0, router, 0)]
        eoc = []
        ic = []
        total_out_ports = 0
        for rid, parameters in reaction_set.iteritems():
            output_port_amount = max(parameters['routing_table'].values()) + 1
            total_out_ports = max(output_port_amount, total_out_ports)
            reaction = self.coder.write_atomic_model('reaction',
                                                     rid,
                                                     [],  # TODO: put correct parameters
                                                     output_port_amount,
                                                     1,
                                                     'Product',
                                                     'Reactant')
            ic.append((router, 0, reaction, 0))
            eoc += [(reaction, port_number, port_number)
                    for port_number in range(output_port_amount)]
            submodels.append(reaction)

        ports = [(0, 'Reactant', 'in')]
        ports += [(port_number, 'Product', 'out') for port_number in range(total_out_ports)]

        return self.coder.write_coupled_model(reaction_set_id,
                                              submodels,
                                              ports,
                                              eic,
                                              eoc,
                                              ic)

    def generate_organelle(self, compartment):

        reaction_sets = self.generate_reaction_sets(compartment)
        print reaction_sets
        space = self.coder.write_atomic_model('space',
                                              compartment.id,
                                              [],
                                              len(reaction_sets),
                                              1,
                                              'Reactants',
                                              'Products')
        submodels = [space] + [model for (model, _, _) in reaction_sets]

        ic = [(space, i, reaction_sets[i][0], 0) for i in range(len(reaction_sets))]
        ic += [(reaction_set, 0, space, 0) for (reaction_set, _, _) in reaction_sets]

        # All membranes use output port 0 to send messages to the space, this is why the
        # coupled output port amount is one less than the max output port amount of its
        # submodels.
        out_port_amount = max([out_port_amount for (_, _, out_port_amount) in reaction_sets])
        ports = [(0, 'Reactants', 'in')]
        ports += [(port_number, 'Product', 'out') for port_number in range(out_port_amount - 1)]

        # if the output port amount is equal to 1, then the model only sends messages to the space
        # and there is no external communication, thus, it must not be linked in the EOC.
        # Because the port 0 of each reaction set goes to the space, the output ports range is
        # [1, .. , output port amount)
        eoc = [(reaction_set, port_number, port_number - 1)
               for (reaction_set, _, output_port_amount) in reaction_sets
               for port_number in range(1, output_port_amount)]

        # if the output port amount is equal to 1, then the model only sends messages to the space
        # and there is no external communication, thus, it must not be linked in the EIC
        eic = [(0, reaction_set, 0)
               for (reaction_set, _, output_port_amount) in reaction_sets
               if output_port_amount > 1]

        return self.coder.write_coupled_model(compartment.id,
                                              submodels,
                                              ports,
                                              eic,
                                              eoc,
                                              ic)
