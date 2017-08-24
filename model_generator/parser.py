#!/usr/bin/python
# -*- coding: utf-8 -*-

# NOTES AND LINKS:
# python SBML parser: https://github.com/linsalrob/PyFBA/blob/master/PyFBA/parse/SBML.py

from bs4 import BeautifulSoup
from collections import defaultdict
import json
import re
# TODO: make all the method camel case


class IlegalCompartmentCombination(Exception):
    pass


class SBMLParser:
    '''
    Parses and retrieve SBML information to generate a CADMIUM P-DEVS model

    '''

    def __init__(self, sbml_file, extra_cellular_id, periplasm_id, cytoplasm_id):
        '''
        SBMLParser constructor

        :param sbml_file: The SBML xml file path.
        :type sbml_file: string
        :return: None
        :rtype: None
        '''

        # Special values
        self.extra_cellular_id = extra_cellular_id
        self.periplasm_id = periplasm_id
        self.cytoplasm_id = cytoplasm_id

        self.unnamed_enzymes_amount = 0

        self.compartments_species = {}
        self.compartments = {}
        self.reactions = {}
        self.enzymes = {}
        self.reaction_locations = {}

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

    def getEnzymesInformation(self):
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
                        'handled_reacions': {rid: self.getReactionParameter(rid)}
                    }
                else:
                    self.enzymes[ehid]['handled_reacions'][rid] = self.getReactionParameter(rid)

        return self.enzymes

    def getReactionIds(self):
        '''
        :param:
        :return: All the SBML reaction IDs without including the biomass reaction
        :rtype: List of SBML reaction IDs
        '''

        return [reaction.get('id')
                for reaction in self.model.findAll('reaction')
                if not self.is_biomass(reaction.get('id'))]

    def getReactionParameter(self, rid):
        '''
        :return: The atomic model parameters ready to be instanciated for the reaction with id rid.
        :rtype: Dictionary of reaction id and parameters
        '''
        print 'getReactionParameter: ', rid

        if rid in self.reactions.keys():
            return self.reactions[rid]

        reaction = self.model.find('reaction', {'id': rid})
        self.reactions[rid] = {
            # TODO(Routing): location is currently not used in the atomic model,
            # but should be used to determine which port to send the metabolites
            'location': self.getReactionLocation(rid),
            'reversible': False if reaction.get('reversible') == 'false' else True,
            'substrate_sctry': self.getReactionSctry(rid, 'listOfReactants'),
            'products_sctry': self.getReactionSctry(rid, 'listOfProducts'),
            'routing_table': self.getReactionRoutingTable(rid),
            'konSTP': self.konSTPs[rid],
            'konPTS': self.konPTSs[rid],
            'koffSTP': self.koffSTPs[rid],
            'koffPTS': self.koffPTSs[rid]
        }

        return self.reactions[rid]

    def getReactionRoutingTable(self, rid):
        species = [s.get('species') for s in self.model.find('reaction', {'id': rid}).findAll('speciesReference')]
        compartment = self.getReactionLocation(rid)['compartment']
        return {s: 0 if compartment == self.specieCompartment(s) else 1 for s in species}

    def getReactionSctry(self, rid, list_name):
        '''
        Calculates stoichiometry numbers of the reaction reactants or products respectively.

        :param: rid, The reaction id to gete the stoichiometry.
        :ptype: String
        :param: list_name, string - One of listOfReactants or listOfProduct in order to retrieve the
        :ptype: String
        :return: The products/reactants stoichiometry
        :rtype: dictionary
        '''

        sctries = self.model.find('reaction', {'id': rid}).find(list_name)

        if sctries is not None:
            return {s.get('species'): 1.0
                    if s.get('stoichiometry') is None
                    else float(s.get('stoichiometry'))
                    for s in sctries.findAll('speciesReference')}
        else:
            return {}

    def newLocation(self, compartment_id, reaction_set):
        return {'compartment': compartment_id, 'reaction_set': reaction_set}

    def getReactionSetIds(self, compartment_id, reaction_set):
        location = self.newLocation(compartment_id, reaction_set)
        reaction_ids = [r.get('id') for r in self.model.findAll('reaction')]
        return [r for r in reaction_ids if self.getReactionLocation(r) == location]

    def getReactionLocation(self, rid):
        # TODO: Get a more general and flexible implementation for this in order to
        # accept different structures. An option is to load a map of str(set<compartment>) as the key
        # and the location as the value. The map can be loaded from a file allowing different
        # compartment combination interpretations.

        if rid in self.reaction_locations.keys():
            return self.reaction_locations[rid]

        reactants = self.getReactionSctry(rid, 'listOfReactants').keys()
        products = self.getReactionSctry(rid, 'listOfProducts').keys()

        compartments = set([self.specieCompartment(s) for s in set(products + reactants)])

        location = None
        if len(compartments) == 1:
            location = self.newLocation(compartments.pop(), 'bulk')

        elif len(compartments) == 2:
            if self.periplasm_id in compartments:  # Periplasm inner or outer membrane

                compartments.discard(self.periplasm_id)
                reaction_sets = {self.extra_cellular_id: 'outer', self.cytoplasm_id: 'inner'}
                location = self.newLocation(self.periplasm_id, reaction_sets[compartments.pop()])
            elif self.cytoplasm_id in compartments:

                compartments.discard(self.cytoplasm_id)
                if self.extra_cellular_id in compartments:  # Periplasm trans membrane
                    location = self.newLocation(self.periplasm_id, 'trans')
                else:  # Organelle membrane to cytoplasm
                    location = self.newLocation(compartments.pop(), 'membrane')

        elif len(compartments) == 3:  # Periplasm trans membrane
            if set([self.extra_cellular_id, self.periplasm_id, self.cytoplasm_id]) == compartments:
                location = self.newLocation(self.periplasm_id, 'trans')

        if location is not None:
            self.reaction_locations[rid] = location
            return location
        else:
            raise IlegalCompartmentCombination('Ilegal compartment combination')

    def specieCompartment(self, specie_id):
        return self.model.find('species', {'id': specie_id}).get('compartment')
