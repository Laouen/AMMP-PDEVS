#!/usr/bin/python
# -*- coding: utf-8 -*-

from SBMLParser import SBMLParser

class ModelStructure:

    def __init__(self, comp_id, internal_reaction_sets, parser, external_reaction_sets=[]):
        self.id = comp_id
        self.space = {}
        self.routing_table = {}
        self.reaction_sets = {}

        structure = ModelStructure(comp_id)

        # Building routing table, each reaction can be reached using the port
        # that the routing table indicates. Each port communicates the space with
        # a different reaction_set
        routing_sets = internal_reaction_sets + external_reaction_sets
        rs_len = range(len(routing_sets))
        structure.routing_table = {(comp_id, r): p for r, p in zip(routing_sets, rs_len)}

        # Building reaction sets
        for reaction_set in internal_reaction_sets:
            reaction_ids = parser.getReactionSetIds(comp_id, reaction_set)
            structure.reaction_sets = {rid: parser.reactions[rid] for rid in reaction_ids}

        # Building space
        reaction_locations = parser.getRelatedReactions(comp_id)
        structure.space = {
            'species': parser.getSpecieByCompartments()[comp_id].keys(),
            'reaction_locations': reaction_locations,
            'enzymes': parser.getRelatedEnzymes(set(reaction_locations.keys()))
        }

class ModelGenerator:
    """
    Generates a model structure from a SBML file using the SBMLParser to get the information
    """

    def __init__(self, sbml_file, extra_cellular_id, periplasm_id, cytoplasm_id, reactions=None, enzymes=None):

        self.parser = SBMLParser(sbml_file, extra_cellular_id, periplasm_id, cytoplasm_id, reactions, enzymes)

        periplasm_sets = ['outer', 'inner', 'trans', 'bulk']
        special_comp_ids = [extra_cellular_id, periplasm_id, cytoplasm_id]

        # TODO: correctly set the cytoplasm and extra cellular space external reaction sets
        self.periplasm = ModelStructure(periplasm_id, periplasm_sets, self.parser)
        self.extra_cellular = ModelStructure(extra_cellular_id, ['bulk'], self.parser)
        self.cytoplasm = ModelStructure(cytoplasm_id, ['bulk'], self.parser)
        self.organelles = [ModelStructure(comp_id, ['membrane', 'bulk'], self.parser)
                           for comp_id in self.parser.getCompartments()
                           if comp_id not in special_comp_ids]
