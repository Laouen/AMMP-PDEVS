#!/usr/bin/python
# -*- coding: utf-8 -*-

## @package pyexample
# Documentation for SBMLParser module.
# This module is a wrapper for lxml: http://lxml.de/xpathxslt.html#the-xpath-method

from lxml import etree
from collections import defaultdict
import json
import re
import gflags
import sys


class SBMLParser:
    ## This class parses, retrieve SBML information and generate a
    #  CADMIUM P-DEVS model

    def __init__(self, sbml_file, sbml_namespace, w3_namespace, biomassIDs=[]):
        ## SBMLParser constructor
        # @param sbml_file: string - The SBML xml file path.
        # @param sbml_namespace: string - The xmlns attribute  of the <sbml> tag.
        # @param w3_namespace: string - The xmlns attribute  of the <body> tags.

        self.sbml_namespace = {'sbml': sbml_namespace}
        self.w3_namespace = {'w3': w3_namespace}
        self.tree = etree.parse(sbml_file)
        self.unnamed_enzymes = 0
        self.biomassIDs = biomassIDs
        self.enzyme_amounts = defaultdict(lambda: 0)
        self.konSTPs = defaultdict(lambda: 0)
        self.konPTSs = defaultdict(lambda: 0)
        self.koffSTPs = defaultdict(lambda: 0)
        self.koffPTSs = defaultdict(lambda: 0)
        self.speciesByCompartments = {}
        self.compartments = {}
        self.reactions = {}
        self.enzymes = {}

    def getNodes(self, tag_name):
        ## Returns all the nodes with tag name = tag_name.
        # @param tag_name: string - The tags name to return.
        # @return [nodes]:

        return self.tree.xpath('//sbml:' + tag_name, namespaces=self.sbml_namespace)

    def getCompartments(self):
        ## Returns a dictionary with the compartment ids as keys and the
        # compartment names as values
        # @return {compartment_id: compartment_name}:

        if not self.compartments == {}:
            return self.compartments

        for compartment in self.getNodes('compartment'):
            self.compartments[compartment.get('id')] = compartment.get('name')

    def getSpecieByCompartments(self):
        ## Returns a dictionary with the compartment ids as keys and dictionaries
        # of the compartment species id and name as value.
        # @return {compartment_id: {specie_id: specie_name}}:

        if self.speciesByCompartments:
            return self.speciesByCompartments

        self.speciesByCompartments = {k: {} for k in self.getCompartments().keys()}

        for specie in self.getNodes('species'):
            id = specie.get('id')
            compartment = specie.get('compartment')
            name = specie.get('name')

            self.speciesByCompartments[compartment][id] = name
        return self.speciesByCompartments

    def getGeneAssociation(self, reaction):
        ## Gets the GENE_ASSOCIATION: text informations from a xml reaction node.
        #  If there is no associated enzymes, it returns an empty list, if there is some associated
        #  enzymes, it returns a list of all of them.
        #  @param reaction: lxml Element - The reaction node from the SBML xml file

        return [x.text.replace('GENE_ASSOCIATION: ', '')
                for x in reaction.xpath('.//w3:p', namespaces=self.w3_namespace)
                if 'GENE_ASSOCIATION' in x.text and any(char.isdigit() for char in x.text)]

    def getEnzymesHandlerIDs(self, reaction):
        ## Returns all the reaction associated enzymes, all the associated enzymes are the enzymes responsables to
        #  handle the reaction.
        #  @param reaction: lxml Element - The reaction node from the SBML xml file

        gene_association = self.getGeneAssociation(reaction)

        if len(gene_association) > 0:
            return self.getEnzymesFromGaneAssociation(gene_association)

        res = ["unnamed_" + str(self.unnamed_enzymes)]
        self.unnamed_enzymes += 1
        return res

    def getEnzymesFromGaneAssociation(self, gene_association):
        ## Returns a list of enzymes from a gene_association logic expresion as the one specified in the SBML files
        #  @param gene_association: string - The gene association logical expresion

        gene_association = map(lambda x: x.replace('(', '[').replace(')', ']').replace(' ', ','), gene_association)
        gene_association = map(lambda x: re.sub(r'([a-z]+[0-9]+)', r'"\1"', x), gene_association)
        gene_association = map(lambda x: re.sub(r'(and|or)', r'"\1"', x), gene_association)
        gene_association = map(lambda x: json.loads(x), gene_association)
        enzymes = map(lambda x: self.parseLogicalGeneAssociation(x), gene_association)
        return [handler for handler_list in enzymes for handler in handler_list]

    def parseLogicalGeneAssociation(self, enzymes):
        ## Parses a gene association expresion of a reaction to get all the responzable enzymes for the reaction
        # @param enzymes: [string] - A list of gene association expresions to parse
        # @return [gene_names]

        if type(enzymes) is unicode or type(enzymes) is str:
            if enzymes in ['and', 'or']:
                return enzymes
            else:
                return [enzymes]

        enzymes = [self.parseLogicalGeneAssociation(expresion) for expresion in enzymes]

        while len(enzymes) >= 3:
            left = enzymes[0]
            operator = enzymes[1]
            right = enzymes[2]
            del enzymes[:3]

            if operator == 'or':
                enzymes = [left + right] + enzymes
            elif operator == 'and':
                product = []
                for i in left:
                    for j in right:
                        product.append(i + "-" + j)
                enzymes = [product] + enzymes

        return enzymes[0]

    def getEnzymes(self):
        ## Returns a map of all the needed enzyme information to intianciate the space models.

        if not self.enzymes == {}:
            return self.enzymes

        if self.reactions == {}:
            self.getReactions()

        for reaction in self.getNodes('reaction'):
            reaction_id = reaction.get('id')
            if reaction_id in self.biomassIDs:
                continue

            for enzyme_handler in self.getEnzymesHandlerIDs(reaction):
                if enzyme_handler not in self.enzymes.keys():
                    self.enzymes[enzyme_handler] = {
                        'id': enzyme_handler,
                        'amount': self.enzyme_amounts[enzyme_handler],
                        'handled_reacions': {reaction_id: self.reactions[reaction_id]}
                    }
                else:
                    self.enzymes[enzyme_handler]['handled_reacions'][reaction_id] = self.reactions[reaction_id]

    def getReactionIDs(self, biomassIDs):
        ## Returns a list with all the SBML reaction IDs
        #  @param biomassIDs: [string] - A list of biomass reaction IDs to not be included in the returned list.
        return [reaction.get('id') for reaction in self.getNodes('reaction')
                if reaction.get('id') not in biomassIDs]

    def getReactions(self):
        ## Returns a map with all the reaction atomic model parameters aready to instanciate them.

        if not self.reactions == {}:
            return self.reactions

        for reaction in self.getNodes('reaction'):
            reaction_id = reaction.get('id')
            if reaction_id in self.biomassIDs:
                continue

            self.reactions[reaction_id] = {
                'reversible': False if reaction.get('reversible') == 'false' else True,
                'substrate_sctry': self.getReactionSctry(reaction, 'listOfReactants'),
                'products_sctry': self.getReactionSctry(reaction, 'listOfProducts'),
                'location': '',  # TODO: implement the method to get the location
                'konSTP': self.konSTPs[reaction_id],
                'konPTS': self.konPTSs[reaction_id],
                'koffSTP': self.koffSTPs[reaction_id],
                'koffPTS': self.koffPTSs[reaction_id]
            }

    def getReactionSctry(self, reaction, list_name):
        ## returns the stoichiometry number of the reaction.
        #  @param reaction: lxml Element - The reaction node to gete the stoichiometry.
        #  @param list_name: string - One of listOfReactants or listOfProduct in order to retrieve the
        #  reactant or product stoichiometry respectively.
        sctries = reaction.xpath('.//sbml:{}/sbml:speciesReference'.format(list_name), namespaces=self.sbml_namespace)

        return {s.get('species'): 1.0
                if s.get('stoichiometry') is None else float(s.get('stoichiometry'))
                for s in sctries}


if __name__ == '__main__':

    gflags.DEFINE_string('sbml_file', None, 'The SBML file path to parse', short_name='f')
    gflags.DEFINE_string('sbml_xmlns', 'http://www.sbml.org/sbml/level2', 'The SBML xmlns attribute value')
    gflags.DEFINE_string('w3_xmlns', 'http://www.sbml.org/sbml/level2', 'The notes body xmlns attribute value')

    gflags.MarkFlagAsRequired('sbml_file')
    FLAGS = gflags.FLAGS

    try:
        argv = FLAGS(sys.argv)  # parse flags
    except gflags.FlagsError, e:
        print '%s\nUsage: %s ARGS\n%s' % (e, sys.argv[0], FLAGS)
        sys.exit(1)

    my_parser = SBMLParser(FLAGS.sbml_file, FLAGS.xmlns)
