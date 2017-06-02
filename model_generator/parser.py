#!/usr/bin/python
# -*- coding: utf-8 -*-
"""@package docstring
Documentation for SBMLParser module.
"""

from lxml import etree
import gflags
import sys

class SBMLParser:
  """ This class parses and retrieve the SBML information to generate a 
      CADMIUM P-DEVS model """

  def __init__(self, sbml_file, namespaces):
    """ SBMLParser constructor

    It parses the sbml_file and saves the namespaces to future use.
    """
    self.namespace = {'sbml': namespaces}
    self.tree = etree.parse(sbml_file)

  def getNodes(self, tag_name):
    """ SBMLParser constructor
    
    It parses the sbml_file and saves the namespaces to future use.
    """
    return self.tree.xpath('//sbml:'+tag_name, namespaces=self.namespace)

if __name__ == '__main__':

  gflags.DEFINE_string('sbml_file', None, 'The SBML file path to parse', short_name='f')
  gflags.DEFINE_string('xmlns', 'http://www.sbml.org/sbml/level2', 'The SBML xmlns attribute value', short_name='x')

  gflags.MarkFlagAsRequired('sbml_file') 
  FLAGS = gflags.FLAGS

  try:
    argv = FLAGS(sys.argv)  # parse flags
  except gflags.FlagsError, e:
    print '%s\nUsage: %s ARGS\n%s' % (e, sys.argv[0], FLAGS)
    sys.exit(1)

  my_parser = SBMLParser(FLAGS.sbml_file,FLAGS.xmlns)
  for specie in my_parser.getNodes('species'):
    print(specie.get("name"))