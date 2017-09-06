#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
NOTES AND LINKS:
python SBML parser: https://github.com/linsalrob/PyFBA/blob/master/PyFBA/parse/SBML.py
"""

import json
import progressbar # sudo pip install progressbar
import re
from bs4 import BeautifulSoup # sudo pip install beautifulsoup, lxml
from collections import defaultdict


class IlegalCompartmentCombination(Exception):
    pass

class SBMLParser:
    """
    Parses and retrieve SBML information to generate a CADMIUM P-DEVS model
    """

    def __init__(self, sbml_file, extra_cellular_id, periplasm_id, cytoplasm_id, reactions={}, enzymes={}):
        """
        class SBMLParser:
        SBMLParser constructor

        :param sbml_file: The SBML xml file path.
        :type sbml_file: string
        :return: None
        :rtype: None
        """

        # Special values
        self.extra_cellular_id = extra_cellular_id
        self.periplasm_id = periplasm_id
        self.cytoplasm_id = cytoplasm_id

        self.unnamed_enzymes_amount = 0

        self.compartments_species = {}
        self.compartments = {}
        self.reactions = reactions
        self.enzymes = enzymes
        self.reaction_locations = {}

        # Not in the SBML model parameters
        self.enzyme_amounts = defaultdict(lambda: 0)
        self.konSTPs = defaultdict(lambda: 0)
        self.konPTSs = defaultdict(lambda: 0)
        self.koffSTPs = defaultdict(lambda: 0)
        self.koffPTSs = defaultdict(lambda: 0)

        self.loadSbmlFile(sbml_file)

        if self.reactions == {} and self.enzymes == {}:
            self.parseReactions()

    def loadSbmlFile(self, sbml_file):
        self.model = BeautifulSoup(open(sbml_file, 'r'), 'xml')

        '''
        TODO: Parse de compartments
        TODO: Parse de reactions
        TODO: Parse de enzymes
        '''

    def isBiomass(self, rid):
        return 'biomass' in rid.lower()

    def getCompartments(self):
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

    def getSpecieByCompartments(self):
        '''
        :return: A dictionary of compartment ids as keys and dictionaries of (specie_id, specie_name) as values 
        :rtype: diccionary<int, dictionary<int, string>>
        '''
        if self.compartment_species:
            return self.compartment_species

        self.compartment_species = {k: {} for k in self.getCompartments().keys()}

        for specie in self.model.findAll('species'):
            id = specie.get('id')
            compartment = specie.get('compartment')
            name = specie.get('name')

            self.compartment_species[compartment][id] = name

        return self.compartment_species

    ##### SPACE ######
    def getRelatedReactionsLocations(self, comp_id):
        '''
        :return: All the reactions that uses species from the compartment and its locations
        :rtype: List of rid (Reaction ID)
        '''

        species = set(self.getSpecieByCompartments()[comp_id].keys())
        return {rid: params['location'] for rid, params in self.reactions 
                if len(species.intersection(params['species'])) > 0}

    def getRelatedEnzymes(self, reactions):
        '''
        :param reactions: A set of reactions.
        :type: set<rid>
        :return: The list of eid (enzyme IDs) that handle the reactions passed as parameter
        :rtype: list<eid>
        '''
        return [eid for eid, enzyme in self.enzymes.items() 
                if len(reactions.intersection(enzyme['handled_reacions'])) > 0]

    def getEnzymeIds(self, reaction):
        '''
        Returns all the reaction associated enzymes, all the associated enzymes are the enzymes responsables to
        handle the reaction.

        :param: reaction:
        :ptype: BeautifulSoup reaction node
        :return: the enzyme handler ids
        :rtype: List of string
        '''

        gene_association = self.getGeneAssociation(reaction)

        if len(gene_association) > 0:
            return self.parseGeneAssociation(gene_association)

        res = ["enzyme_" + str(self.unnamed_enzymes_amount)]
        self.unnamed_enzymes_amount += 1
        return res

    def getGeneAssociation(self, reaction):
        '''
        Gets the GENE_ASSOCIATION: text informations from a xml reaction node.
        If there is no associated enzymes, it returns an empty list, if there 
        is some associated enzymes, it returns a list of all of them.

        :param reaction: The reaction node from the SBML xml file
        :type: lxml Element.
        :return: A list of logical gene associations.
        :rtype: List<string>
        '''
        return [x.text.replace('GENE_ASSOCIATION: ', '')
                for x in reaction.notes.body.findAll('p')
                if 'GENE_ASSOCIATION' in x.text and any(char.isdigit() for char in x.text)]

    def parseGeneAssociation(self, gene_association):
        '''
        :param gene_association: The gene association logical expresion
        :type: string
        :return: A list of enzyme ids from a gene_association logic expresion as the one specified in the SBML files
        :rtype: List<string>
        '''
        gene_association = map(lambda x: x.replace('(', '[').replace(')', ']').replace(' ', ','), gene_association)
        gene_association = map(lambda x: re.sub(r'([a-z]+[0-9]+)', r'"\1"', x), gene_association)
        gene_association = map(lambda x: re.sub(r'(and|or)', r'"\1"', x), gene_association)
        gene_association = map(lambda x: json.loads(x), gene_association)
        enzymes = map(lambda x: self.parseLogicalGeneAssociation(x), gene_association)
        return [handler for handler_list in enzymes for handler in handler_list]

    def parseLogicalGeneAssociation(self, enzymes):
        '''
        Parses a gene association expresion of a reaction to get all the responzable enzymes for the reaction

        :param enzymes: A list of gene association expresions to parse
        :type: list<string>
        :return: All the gene names associated to the enzymes
        :rtype: list<string>
        '''

        # base cases
        if type(enzymes) in [unicode, str]:
            if enzymes in ['and', 'or']:
                return enzymes
            else:
                return [enzymes]

        # recursive cases
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
                for l in left:
                    for r in right:
                        product.append(l + '-' + r)
                enzymes = [product] + enzymes

        return enzymes[0]

    def parseEnzymes(self, reaction):
        '''
        Parses all the reaction enzymes and saves their information to be used 
        by spaces retrieved information: eid (Enzyme ID), handled reactions.
        '''

        rid = reaction.get('id')
        if self.isBiomass(rid):
            return

        for eid in self.getEnzymeIds(reaction):
            if eid not in self.enzymes.keys():
                self.enzymes[eid] = {
                    'id': eid,
                    # TODO: the amount should be in the model generator and 
                    # should be per compartments
                    'amount': self.enzyme_amounts[eid],
                    'handled_reactions': [rid]
                }
            else:
                self.enzymes[eid]['handled_reactions'].append(rid)

    ##### REACTIONS ######
    ## Implemented optimized methods where all re resources are acceded at once

    def getReactionSetIds(self, compartment_id, reaction_set):
        location = self.newLocation(compartment_id, reaction_set)
        return [r for r, p in self.reactions.items() if p['location'] == location]

    def parseReactions(self):
        '''
        Parses the reaction informations and construct the reaction parameters 
        for the model. 
        '''
        print '[Parser] Start parsing reactions.'
        reactions = self.model.findAll('reaction')
        
        bar = progressbar.ProgressBar(maxval=len(reactions))
        bar.start()
        done = 0
        for reaction in reactions:
            self.parseReaction(reaction)
            bar.update(done)
            done += 1
        
        print '[Parser] End parsing reactions.'

    def newLocation(self, compartment_id, reaction_set):
        return {'compartment': compartment_id, 'reaction_set': reaction_set}

    def specieCompartment(self, specie_id):
        return self.model.find('species', {'id': specie_id}).get('compartment')

    def getReactionRoutingTable(self, reaction, comp_id):
        species = [s.get('species') for s in reaction.findAll('speciesReference')]
        return {s: 0 if comp_id == self.specieCompartment(s)
                else 1 if self.cytoplasm_id == self.specieCompartment(s)
                else 2 for s in species}

    def parseStoichiometry(self, reaction):
        stoichiometry = {}

        # parsing reactants
        sctries = reaction.find('listOfReactants')
        if sctries is not None:
            reactants = {
                s.get('species'): 1.0
                if s.get('stoichiometry') is None
                else float(s.get('stoichiometry'))
                for s in sctries.findAll('speciesReference')
            }
        else:
            reactants = {}

        # parsing products
        sctries = reaction.find('listOfProducts')
        if sctries is not None:
            products = {
                s.get('species'): 1.0
                if s.get('stoichiometry') is None
                else float(s.get('stoichiometry'))
                for s in sctries.findAll('speciesReference')
            }
        else:
            products = {}
        
        stoichiometry['listOfReactants'] = reactants
        stoichiometry['listOfProducts'] = products
        return stoichiometry

    def getLocation(self, species):
        '''
        TODO: Get a more general and flexible implementation for this in order to
        accept different structures. An option is to load a map of str(set<compartment>) as the key
        and the location as the value. The map can be loaded from a file allowing different
        compartment combination interpretations.
        '''

        compartments = set([self.specieCompartment(s) for s in set(species)])

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

        if location is None:
            raise IlegalCompartmentCombination('Ilegal compartment combination ' + str(compartments))
        
        return location

    def parseReaction(self, reaction):
        rid = reaction.get('id')
        if self.isBiomass(rid):
            return

        stoichiometry = self.parseStoichiometry(reaction)
        species = stoichiometry['listOfReactants'].keys() + stoichiometry['listOfProducts'].keys()
        location = self.getLocation(species)
        routing_table = self.getReactionRoutingTable(reaction,
                                                     location['compartment'])
        parameters = {
            # TODO(Routing): location is currently not used in the atomic model,
            # but should be used to determine which port to send the metabolites
            'location': location,
            'reversible': False if reaction.get('reversible') == 'false' else True,
            'species': species,
            'substrate_sctry': stoichiometry['listOfReactants'],
            'products_sctry': stoichiometry['listOfProducts'],
            'routing_table': routing_table,
            'konSTP': self.konSTPs[rid],
            'konPTS': self.konPTSs[rid],
            'koffSTP': self.koffSTPs[rid],
            'koffPTS': self.koffPTSs[rid]
        }

        self.parseEnzymes(reaction)
        self.reactions[rid] = parameters