#!/usr/bin/python
# -*- coding: utf-8 -*-

import gflags
import sys
import pickle

from ModelGenerator import ModelGenerator

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

    reactions = None
    enzymes = None

    '''
    with open("pickles/reactions.pickle", "r") as reaction_file:
        reactions = pickle.load(reaction_file)

    with open("pickles/enzymes.pickle", "r") as enzymes_file:
        enzymes = pickle.load(enzymes_file)
    
    '''

    generator = ModelGenerator(FLAGS.sbml_file,
                               FLAGS.extra_cellular,
                               FLAGS.periplasm,
                               FLAGS.cytoplasm,
                               reactions=reactions,
                               enzymes=enzymes)
    """
    with open("pickles/reactions.pickle", "wb") as reaction_file:
        pickle.dump(generator.parser.reactions, reaction_file)

    with open("pickles/enzymes.pickle", "wb") as enzymes_file:
        pickle.dump(generator.parser.enzymes, enzymes_file)
    """
    generator.generate_top()
