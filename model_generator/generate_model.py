#!/usr/bin/python
# -*- coding: utf-8 -*-

import gflags
import sys
from modelGenerator import ModelGenerator

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
