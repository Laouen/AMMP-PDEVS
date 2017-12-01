#!/usr/bin/python
# -*- coding: utf-8 -*-

from SBMLTestGenerator import SBMLTestGenerator
import gflags  # sudo pip install gflags
import sys
from tqdm import tqdm


def main(FLAGS):

    models_range = len(FLAGS.reaction_amount)
    assert models_range == len(FLAGS.species_amount) and models_range == len(FLAGS.max_stoichimetry_elements)

    for i in tqdm(range(models_range)):
        reaction_amount = int(FLAGS.reaction_amount[i])
        species_amount = int(FLAGS.species_amount[i])
        max_stoichimetry_elements = int(FLAGS.max_stoichimetry_elements[i])

        reaction_amounts = {
            ('c', 'bulk'): reaction_amount,
            ('e', 'bulk'): reaction_amount,
            ('p', 'bulk'): reaction_amount,
            ('p', 'inner'): reaction_amount,
            ('p', 'outer'): reaction_amount,
            ('p', 'trans'): reaction_amount
        }

        compartments = ['c', 'e', 'p']
        species_amounts = {cid: species_amount for cid in compartments}

        model_id = '_'.join([
            FLAGS.model_name,
            str(reaction_amount),
            str(species_amount),
            str(max_stoichimetry_elements)
        ])

        generator = SBMLTestGenerator(model_id + '.xml',
                                      model_id,
                                      compartments,
                                      reaction_amounts,
                                      species_amounts,
                                      max_stoichimetry_elements)
        generator.create_listOfCompartments()
        generator.create_listOfSpecies()
        generator.create_listOfReactions()
        generator.save_xml()


if __name__ == '__main__':

    gflags.DEFINE_list('reaction_amount', [10], 'A list with the amounts of reactions in each compartment',
                       short_name='r')
    gflags.DEFINE_list('species_amount', [100], 'A list with the amounts of species in each compartment',
                       short_name='s')
    gflags.DEFINE_list('max_stoichimetry_elements', [10], 'The maximum amount of Product and '
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
