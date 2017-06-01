#!/usr/bin/python
# -*- coding: utf-8 -*-

from lxml import etree
import gflags
import sys

if __name__ == '__main__':

  gflags.DEFINE_string('smbl_file', None, 'The SBML file path to parse', short_name='f')

  gflags.MarkFlagAsRequired('smbl_file') 
  FLAGS = gflags.FLAGS

  try:
    argv = FLAGS(sys.argv)  # parse flags
  except gflags.FlagsError, e:
    print '%s\nUsage: %s ARGS\n%s' % (e, sys.argv[0], FLAGS)
    sys.exit(1)

  tree = etree.parse(FLAGS.smbl_file)
  print(tree.xpath("/model"))
  for specie in tree.xpath("/model/listOfSpecies/species"):
    print(specie.id)