#!/usr/bin/python
# -*- coding: utf-8 -*-

from SBMLTestGenerator import SBMLTestGenerator
import gflags  # sudo pip install gflags
import sys
from tqdm import tqdm


def main(FLAGS):

    for i in tqdm(range(len(FLAGS.reaction_amount))):

        reaction_amount = int(FLAGS.reaction_amount[i])

        reaction_amounts = {
            ('c', 'bulk'): reaction_amount,
            ('e', 'bulk'): reaction_amount,
            ('p', 'bulk'): reaction_amount,
            ('p', 'inner'): reaction_amount,
            ('p', 'outer'): reaction_amount,
            ('p', 'trans'): reaction_amount
        }

        compartments = ['c', 'e', 'p']

        model_id = '_'.join([
            FLAGS.model_name,
            "%03d" % reaction_amount,  # format number with three digits, 1 -> 001, 10 -> 010
            "%03d" % FLAGS.stoichimetry_elements
        ])

        generator = SBMLTestGenerator(model_id + '.xml',
                                      model_id,
                                      compartments,
                                      reaction_amounts,
                                      FLAGS.stoichimetry_elements)
        generator.save_xml()


if __name__ == '__main__':

    gflags.DEFINE_list('reaction_amount', [1], 'A list with the amounts of reactions in each compartment', short_name='r')
    gflags.DEFINE_integer('stoichimetry_elements', 1, 'The number of Product and Reactant metabolites a reaction can have', short_name='m')
    gflags.DEFINE_string('model_name', 'process_time_test_', 'The SBML model name', short_name='n')

    FLAGS = gflags.FLAGS

    try:
        argv = FLAGS(sys.argv)  # parse flags
    except gflags.FlagsError as e:
        print('%s\nUsage: %s ARGS\n%s' % (e, sys.argv[0], FLAGS))
        sys.exit(1)

    main(FLAGS)
