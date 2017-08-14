#!/usr/bin/python
# -*- coding: utf-8 -*-

## @package pyexample
# Documentation for SBMLParser module.
# This module is a wrapper for lxml: http://lxml.de/xpathxslt.html#the-xpath-method


# NOTES AND LINKS:
# python SBML parser: https://github.com/linsalrob/PyFBA/blob/master/PyFBA/parse/SBML.py

from bs4 import BeautifulSoup
from collections import defaultdict
import json
import re
import gflags
import sys


class SBMLParser:
    '''
    Parses and retrieve SBML information to generate a CADMIUM P-DEVS model

    '''

    def __init__(self, sbml_file):
        '''
        SBMLParser constructor

        :param sbml_file: The SBML xml file path.
        :type sbml_file: string
        :return: None
        :rtype: None
        '''

        self.unnamed_enzymes_amount = 0

        self.compartments_species = {}
        self.compartments = {}
        self.reactions = {}
        self.enzymes = {}

        # Not in the SBML model parameters
        self.enzyme_amounts = defaultdict(lambda: 0)
        self.konSTPs = defaultdict(lambda: 0)
        self.konPTSs = defaultdict(lambda: 0)
        self.koffSTPs = defaultdict(lambda: 0)
        self.koffPTSs = defaultdict(lambda: 0)

        self.parse_sbml_file(sbml_file)

    def parse_sbml_file(self, sbml_file):
        self.model = BeautifulSoup(open(sbml_file, 'r'), 'xml')

        # Parse de compartments
        # Parse de reactions
        # Parse de enzymes

    def is_biomass(eid):
        return 'biomass' in eid.lower()

    def get_compartments(self):
        '''
        Returns a dictionary with the compartment ids as keys and the
        compartment names as values
        :return: compartment_id
        :rtype: dictionary
        '''

        if not self.compartments == {}:
            return self.compartments

        for compartment in self.model.findAll('compartment'):
            self.compartments[compartment.get('id')] = compartment.get('name')

    def get_specie_by_compartments(self):
        ## Returns a dictionary with the compartment ids as keys and dictionaries
        # of the compartment species id and name as value.
        # @return {compartment_id: {specie_id: specie_name}}:

        if self.compartment_species:
            return self.compartment_species

        self.compartment_species = {k: {} for k in self.getCompartments().keys()}

        for specie in self.model.findAll('species'):
            id = specie.get('id')
            compartment = specie.get('compartment')
            name = specie.get('name')

            self.compartment_species[compartment][id] = name

        return self.compartment_species

    def get_gene_association(self, reaction):
        ## Gets the GENE_ASSOCIATION: text informations from a xml reaction node.
        #  If there is no associated enzymes, it returns an empty list, if there is some associated
        #  enzymes, it returns a list of all of them.
        #  @param reaction: lxml Element - The reaction node from the SBML xml file
        #  @return: A list of logical gene associations

        return [x.text.replace('GENE_ASSOCIATION: ', '')
                for x in reaction.notes.body.findAll('p')
                if 'GENE_ASSOCIATION' in x.text and any(char.isdigit() for char in x.text)]

    def get_enzymes_handler_ids(self, reaction):
        '''
        Returns all the reaction associated enzymes, all the associated enzymes are the enzymes responsables to
        handle the reaction.
        :param: reaction:
        :ptype: BeautifulSoup reaction node
        :return: the enzyme handler ids
        :rtype: List of string
        '''

        gene_association = self.get_gene_association(reaction)

        if len(gene_association) > 0:
            return self.get_enzymes_from_gane_association(gene_association)

        res = ["unnamed_" + str(self.unnamed_enzymes_amount)]
        self.unnamed_enzymes_amount += 1
        return res

    def get_enzymes_from_gane_association(self, gene_association):
        ## Returns a list of enzymes from a gene_association logic expresion as the one specified in the SBML files
        #  @param gene_association: string - The gene association logical expresion

        gene_association = map(lambda x: x.replace('(', '[').replace(')', ']').replace(' ', ','), gene_association)
        gene_association = map(lambda x: re.sub(r'([a-z]+[0-9]+)', r'"\1"', x), gene_association)
        gene_association = map(lambda x: re.sub(r'(and|or)', r'"\1"', x), gene_association)
        gene_association = map(lambda x: json.loads(x), gene_association)
        enzymes = map(lambda x: self.parse_logical_gene_association(x), gene_association)
        return [handler for handler_list in enzymes for handler in handler_list]

    def parse_logical_gene_association(self, enzymes):
        ## Parses a gene association expresion of a reaction to get all the responzable enzymes for the reaction
        # @param enzymes: [string] - A list of gene association expresions to parse
        # @return [gene_names]

        # base cases
        if type(enzymes) in [unicode, str]:
            if enzymes in ['and', 'or']:
                return enzymes
            else:
                return [enzymes]

        # recursive cases
        enzymes = [self.parse_logical_gene_association(expresion) for expresion in enzymes]

        while len(enzymes) >= 3:
            left = enzymes[0]
            operator = enzymes[1]
            right = enzymes[2]
            del enzymes[:3]

            if operator == 'or':
                enzymes = [left + right] + enzymes
            elif operator == 'and':
                product = []
                for l in left:
                    for r in right:
                        product.append(l + '-' + r)
                enzymes = [product] + enzymes

        return enzymes[0]

    def get_enzymes(self):
        ## Returns a map of all the needed enzyme information to intianciate the space models.

        if not self.enzymes == {}:
            return self.enzymes

        for reaction in self.model.findAll('reaction'):
            rid = reaction.get('id')
            if self.is_biomass(rid):
                continue

            for ehid in self.get_enzymes_handler_ids(reaction):
                if ehid not in self.enzymes.keys():
                    self.enzymes[ehid] = {
                        'id': ehid,
                        'amount': self.enzyme_amounts[ehid],
                        'handled_reacions': {rid: self.get_reactions()[rid]}
                    }
                else:
                    self.enzymes[ehid]['handled_reacions'][rid] = self.get_reactions()[rid]

    def get_reaction_ids(self):
        ## Returns a list with all the SBML reaction IDs
        #  @param biomassIDs: [string] - A list of biomass reaction IDs to not be included in the returned list.

        return [reaction.get('id')
                for reaction in self.model.findAll('reaction')
                if not self.is_biomass(reaction.get('id'))]

    def get_reactions(self):
        '''
        Returns a map with all the reaction atomic model parameters ready to be instanciated.
        :return: self.reactions
        :rtype: dictionary of reaction id and parameters
        '''
        if not self.reactions == {}:
            return self.reactions

        for reaction in self.model.findAll('reaction'):
            rid = reaction.get('id')
            if self.is_biomass(rid):
                continue

            self.reactions[rid] = {
                'reversible': False if reaction.get('reversible') == 'false' else True,
                'substrate_sctry': self.get_reaction_sctry(rid, 'listOfReactants'),
                'products_sctry': self.get_reaction_sctry(rid, 'listOfProducts'),
                'location': self.get_reaction_location(rid),
                'konSTP': self.konSTPs[rid],
                'konPTS': self.konPTSs[rid],
                'koffSTP': self.koffSTPs[rid],
                'koffPTS': self.koffPTSs[rid]
            }

    # TODO: implement this method
    def get_reaction_location(self, rid):

        return ''

    def get_reaction_sctry(self, rid, list_name):
        '''
        Calculates stoichiometry numbers of the reaction reactants or products respectively.

        :param: rid, The reaction id to gete the stoichiometry.
        :ptype: String
        :param: list_name, string - One of listOfReactants or listOfProduct in order to retrieve the
        :ptype: String
        :return: The products/reactants stoichiometry
        :rtype: dictionary
        '''

        sctries = self.model.sbml.model \
            .find('reaction', {'id': rid}) \
            .find(list_name) \
            .findAll('speciesReference')

        return {s.get('species'): 1.0
                if s.get('stoichiometry') is None
                else float(s.get('stoichiometry'))
                for s in sctries}


if __name__ == '__main__':

    gflags.DEFINE_string('sbml_file', None, 'The SBML file path to parse', short_name='f')

    gflags.MarkFlagAsRequired('sbml_file')
    FLAGS = gflags.FLAGS

    try:
        argv = FLAGS(sys.argv)  # parse flags
    except gflags.FlagsError, e:
        print '%s\nUsage: %s ARGS\n%s' % (e, sys.argv[0], FLAGS)
        sys.exit(1)

    my_parser = SBMLParser(FLAGS.sbml_file)
