#!/usr/bin/python
# -*- coding: utf-8 -*-

from SBMLTestGenerator import SBMLTestGenerator
import gflags  # sudo pip install gflags
import sys


def main(FLAGS):

    reaction_amounts = {
        ('c', 'bulk'): FLAGS.reaction_amount,
        ('e', 'bulk'): FLAGS.reaction_amount,
        ('p', 'bulk'): FLAGS.reaction_amount,
        ('p', 'inner'): FLAGS.reaction_amount,
        ('p', 'outer'): FLAGS.reaction_amount,
        ('p', 'trans'): FLAGS.reaction_amount
    }

    compartments = ['c', 'e', 'p']
    species_amount = {cid: FLAGS.species_amount for cid in compartments}

    generator = SBMLTestGenerator(FLAGS.model_name + '.xml',
                                  FLAGS.model_name,
                                  compartments,
                                  reaction_amounts,
                                  species_amount,
                                  FLAGS.max_stoichimetry_elements)
    generator.create_listOfCompartments()
    generator.create_listOfSpecies()
    generator.create_listOfReactions()
    generator.save_xml()


if __name__ == '__main__':

    gflags.DEFINE_string('reaction_amount', 10, 'The amount of reactions in each compartment', short_name='r')
    gflags.DEFINE_string('species_amount', 100, 'The amount of species in each compartment', short_name='s')
    gflags.DEFINE_string('max_stoichimetry_elements', 10, 'The maximum amount of Product and '
                         'Reactant metabolites a reaction can have', short_name='m')
    gflags.DEFINE_string('model_name', None, 'The SBML model name', short_name='n')

    gflags.MarkFlagAsRequired('model_name')
    FLAGS = gflags.FLAGS

    try:
        argv = FLAGS(sys.argv)  # parse flags
    except gflags.FlagsError, e:
        print '%s\nUsage: %s ARGS\n%s' % (e, sys.argv[0], FLAGS)
        sys.exit(1)

    main(FLAGS)
