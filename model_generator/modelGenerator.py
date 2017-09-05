#!/usr/bin/python
# -*- coding: utf-8 -*-

from parser import SBMLParser

class ModelGenerator:
    '''
    Generates a model structure from a SBML file using the SBMLParser to get the information
    '''

    def __init__(self, sbml_file, extra_cellular_id, periplasm_id, cytoplasm_id):

        self.parser = SBMLParser(sbml_file, extra_cellular_id, periplasm_id, cytoplasm_id)
        self.periplasm = {'id': periplasm_id}
        self.extra_cellular = {'id': extra_cellular_id}
        self.cytoplasm = {'id': cytoplasm_id}
        self.organelles = []

    def generateStructure(self):
        self.generatePeriplasmStructure()

    def generatePeriplasmStructure(self):
        # All the membrane and the bulk are reaction set that sends/receives
        # metabolites to compartment spaces.
        reaction_sets = ['outer', 'inner', 'trans', 'bulk']

        # Building routing table for periplasm, each reaction can be reached using the port
        # that the routing table indicates. Each port communicates the periplasm space with
        # a different membrane
        next_port = 0
        comp_id = self.periplasm['id']
        self.periplasm['routing_table'] = {}
        for reaction_set in reaction_sets:
            self.periplasm['routing_table'][(comp_id, reaction_set)] = next_port
            next_port += 1

        # Building the periplasm reaction_sets structures
        self.periplasm['reaction_sets'] = {}
        for reaction_set in reaction_sets:
            reaction_ids = self.parser.getReactionSetIds(comp_id, reaction_set)
            reaction_set = {rid: self.parser.reactions[rid] for rid in reaction_ids}
            self.periplasm['reaction_sets'][reaction_set] = reaction_set

        # Building the space structure
        reaction_locations = self.parser.getRelatedReactions(comp_id)
        self.periplasm['space'] = {
            'species': self.parser.getSpecieByCompartments()[comp_id].keys(),
            'reaction_locations': reaction_locations,
            'enzymes': self.parser.getRelatedEnzymes(set(reaction_locations.keys()))
        }