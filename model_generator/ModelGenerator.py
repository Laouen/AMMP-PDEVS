#!/usr/bin/python
# -*- coding: utf-8 -*-

from SBMLParser import SBMLParser


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
        print comp_id
        for reaction_set in internal_reaction_sets:
            print reaction_set
            reaction_ids = parser.get_reaction_set_ids(comp_id, reaction_set)
            self.reaction_sets = {rid: parser.reactions[rid] for rid in reaction_ids}

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

    def __init__(self, sbml_file, extra_cellular_id, periplasm_id, cytoplasm_id,
                 reactions=None, enzymes=None):
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
        self.parser = SBMLParser(sbml_file, extra_cellular_id, periplasm_id,
                                 cytoplasm_id, reactions, enzymes)

        periplasm_sets = ['outer', 'inner', 'trans']
        special_comp_ids = [extra_cellular_id, periplasm_id, cytoplasm_id]

        # TODO: correctly set the cytoplasm and extra cellular space external reaction sets
        self.periplasm = ModelStructure(periplasm_id, self.parser, periplasm_sets)
        self.extra_cellular = ModelStructure(extra_cellular_id, self.parser)
        self.cytoplasm = ModelStructure(cytoplasm_id, self.parser)
        self.organelles = [ModelStructure(comp_id, self.parser, ['membrane'])
                           for comp_id in self.parser.get_compartments()
                           if comp_id not in special_comp_ids]
