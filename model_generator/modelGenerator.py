#!/usr/bin/python
# -*- coding: utf-8 -*-

from parser import SBMLParser
import gflags
import sys


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
            self.periplasm['reaction_sets'][reaction_set] = self.getReactionSetStructure(reaction_ids)

        # TODO: implement the space parameters

    def getReactionSetStructure(self, reaction_ids):

        reaction_set = {rid: self.parser.reactions[rid] for rid in reaction_ids}
        return reaction_set


if __name__ == '__main__':

    gflags.DEFINE_string('sbml_file', None, 'The SBML file path to parse', short_name='f')
    gflags.DEFINE_string('extra_cellular', 'e', 'The extra cellular space ID in the SBML file', short_name='e')
    gflags.DEFINE_string('cytoplasm', 'c', 'The cytoplasm space ID in the SBML file', short_name='c')
    gflags.DEFINE_string('periplasm', 'p', 'The periplasm space ID in the SBML file', short_name='p')

    gflags.MarkFlagAsRequired('sbml_file')
    FLAGS = gflags.FLAGS

    try:
        argv = FLAGS(sys.argv)  # parse flags
    except gflags.FlagsError, e:
        print '%s\nUsage: %s ARGS\n%s' % (e, sys.argv[0], FLAGS)
        sys.exit(1)

    generator = ModelGenerator(FLAGS.sbml_file,
                               FLAGS.extra_cellular,
                               FLAGS.periplasm,
                               FLAGS.cytoplasm)

    print 'start geting stoichiometries'
    generator.parser.parseStoichiometries()
    print 'end geting stoichiometries'
    generator.generateStructure()
