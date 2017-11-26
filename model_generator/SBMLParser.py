#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
NOTES AND LINKS:
python SBML parser: https://github.com/linsalrob/PyFBA/blob/master/PyFBA/parse/SBML.py
"""

import json
import progressbar  # sudo pip install progressbar
import re
from bs4 import BeautifulSoup  # sudo pip install beautifulsoup, lxml
from collections import defaultdict
from constants import *
from json import JSONEncoder
import copy


class ReactionParametersDecoder(JSONEncoder):
    def default(self, o):
        print "Reaction Parameter Encoder: ", o
        res = {k: v for k, v in o.iteritems()}
        for key in res.keys():
            res[key]['location'] = dict(res[key]['location'].__dict__)
        return res


class SBMLParserEncoder(JSONEncoder):
    def default(self, o):

        res = {
            'extra_cellular_id': o.extra_cellular_id,
            'periplasm_id': o.periplasm_id,
            'cytoplasm_id': o.cytoplasm_id,
            'unnamed_enzymes_amount': o.unnamed_enzymes_amount,
            'reactions': copy.deepcopy(o.reactions),
            'enzymes': dict(o.enzymes),
        }

        for key in res['reactions'].keys():
            res['reactions'][key]['location'] = dict(res['reactions'][key]['location'].__dict__)

        return res


def not_empty_intersection(a, b):
    return (set(a) & set(b)) != set()


class IllegalCompartmentCombination(Exception):
    pass


class Location:
    """
    Represent a specie location composed by the compartment id (cid attribute) and the
    reaction set name (rsn attribute) where the specie belongs
    """
    def __init__(self, cid, rsn):
        """

        :param cid: The compartment id
        :type cid: str
        :param rsn: Te reaction set name
        :type rsn: str
        """
        self.cid = cid
        self.rsn = rsn

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return other.__dict__ == self.__dict__
        else:
            return False

    def as_tuple(self):
        return tuple([self.cid, self.rsn])


class SBMLParser:
    """
    Parses and retrieve SBML information to generate a CADMIUM P-DEVS model
    """

    def __init__(self,
                 sbml_file=None,
                 extra_cellular_id='e',
                 periplasm_id='p',
                 cytoplasm_id='c',
                 json_model=None):
        """
        class SBMLParser:
        SBMLParser constructor

        :param sbml_file: The SBML xml file path.
        :type sbml_file: string
        :param json_model: The parser exported as json. Optional, used to avoid re parsing
        :type json_mode: dict
        :return: None
        :rtype: None
        """

        self.model = None

        # Special values
        self.extra_cellular_id = extra_cellular_id
        self.periplasm_id = periplasm_id
        self.cytoplasm_id = cytoplasm_id

        self.unnamed_enzymes_amount = 0

        self.compartment_species = None
        self.compartments = None
        self.reactions = None
        self.enzymes = None

        # Not in the SBML model parameters
        # TODO: refactor default dict values to be set in the class init method parameters
        self.enzyme_amounts = defaultdict(lambda: 100)
        self.konSTPs = defaultdict(lambda: 0.8)
        self.konPTSs = defaultdict(lambda: 0.8)
        self.koffSTPs = defaultdict(lambda: 0.8)
        self.koffPTSs = defaultdict(lambda: 0.8)
        self.metabolite_amounts = defaultdict(lambda: 100)
        self.rates = defaultdict(lambda: '0:0:0:1')
        self.reject_times = defaultdict(lambda: '0:0:0:1')
        self.interval_times = defaultdict(lambda: '0:0:0:1')

        self.load_sbml_file(sbml_file)

        if json_model is None:
            self.parse_reactions()
        else:
            self.load_json_model(json_model)

    def load_json_model(self, json_model):
        self.extra_cellular_id = json_model['extra_cellular_id']
        self.periplasm_id = json_model['periplasm_id']
        self.cytoplasm_id = json_model['cytoplasm_id']
        self.unnamed_enzymes_amount = json_model['unnamed_enzymes_amount']
        self.reactions = dict(json_model['reactions'])
        self.enzymes = dict(json_model['enzymes'])

        for key in self.reactions.keys():
            cid = self.reactions[key]['location']['cid']
            rsn = self.reactions[key]['location']['rsn']
            self.reactions[key]['location'] = Location(cid, rsn)

    def load_sbml_file(self, sbml_file):
        """
        Loads the xml sbml file from the file path passed as parameter
        :param sbml_file: The path to the xml file to parse as sbml
        :type sbml_file: str
        """
        if sbml_file is not None:
            self.model = BeautifulSoup(open(sbml_file, 'r'), 'xml')

    def get_compartments(self):
        """

        :return: A dictionary with the cid (compartment IDs) and compartment names
        :rtype: dict[str, str]
        """

        if self.compartments is not None:
            return self.compartments

        self.compartments = {}
        for compartment in self.model.findAll('compartment'):
            self.compartments[compartment.get('id')] = compartment.get('name')

        return self.compartments

    def parse_compartments_species(self):
        """
        :return: A dictionary of compartment ids as keys and dictionaries of
        (specie_id, specie_name) as values
        :rtype: dict[int, dict[int, string]]
        """

        if self.compartment_species is not None:
            return self.compartment_species

        self.compartment_species = {k: {} for k in self.get_compartments().keys()}

        for specie in self.model.findAll('species'):
            specie_id = specie.get('id')
            compartment = specie.get('compartment')
            name = specie.get('name')

            self.compartment_species[compartment][specie_id] = name

        return self.compartment_species

    def get_reaction_parameters(self, comp_id):
        """
        :return: All the reaction parameters for the reactions that uses species from the
        compartment
        :rtype: list[str]
        """
        compartment_species = self.parse_compartments_species()[comp_id].keys()
        return {rid: reaction_parameters
                for rid, reaction_parameters in self.reactions.items()
                if not_empty_intersection(compartment_species, reaction_parameters['species'])}

    def get_enzymes(self, reactions):
        """
        :param reactions: A set of rid (reaction IDs).
        :type: set[str]
        :return: The dictionary of eid (enzyme IDs) as keys and enzyme parameters as values
        of the enzymes that handle the reactions passed as parameter. All enzyme that handles at
        least one of the reactions is added.
        :rtype: dict[str, dict[str, any]]
        """
        return {eid: enzyme_parameters
                for eid, enzyme_parameters in self.enzymes.items()
                if not_empty_intersection(enzyme_parameters['handled_reactions'], reactions)}

    def get_enzyme_ids(self, reaction):
        """
        Returns all the reaction associated enzymes, all the associated enzymes are the enzymes
        responsables to handle the reaction.

        :param: reaction:
        :ptype: bs4.element.Tag
        :return: the enzyme handler ids
        :rtype: list[str]
        """

        gene_association = self.extract_gene_association(reaction)

        if len(gene_association) > 0:
            return self.parse_gene_association(gene_association)

        res = ["enzyme_" + str(self.unnamed_enzymes_amount)]
        self.unnamed_enzymes_amount += 1
        return res

    def parse_gene_association(self, gene_association):
        """
        :param gene_association: The gene association logical expresion
        :type: string
        :return: A list of enzyme ids from a gene_association logic expresion as the one specified
        in the SBML files
        :rtype: List[str]
        """
        gene_association = map(lambda x: x.replace('(', '[').replace(')', ']').replace(' ', ','),
                               gene_association)
        gene_association = map(lambda x: re.sub(r'([a-z]+[0-9]+)', r'"\1"', x),
                               gene_association)
        gene_association = map(lambda x: re.sub(r'(and|or)', r'"\1"', x),
                               gene_association)
        gene_association = map(lambda x: json.loads(x),
                               gene_association)

        # all_enzymes is a list of enzyme lists and the return statement is the one who flattens
        # the all_enzymes in order to correctly return an enzyme list
        all_enzymes = map(lambda x: self.parse_logical_gene_association(x), gene_association)
        return [enzyme for enzyme_list in all_enzymes for enzyme in enzyme_list]

    def parse_logical_gene_association(self, enzymes):
        """
        Parses a gene association expresion of a reaction to get all the responsible
        enzymes for the reaction

        :param enzymes: A list of gene association expresions to parse
        :type: list[str]
        :return: All the gene names associated to the enzymes
        :rtype: list[str]
        """

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

    def parse_enzymes(self, reaction):
        """
        Parses all the reaction enzymes and saves their information to be used
        by spaces retrieved information: eid (Enzyme ID), handled reactions.
        """

        rid = reaction.get('id')
        if self.is_biomass(rid) is True:
            return

        for eid in self.get_enzyme_ids(reaction):
            if eid not in self.enzymes.keys():
                self.enzymes[eid] = {
                    'id': eid,
                    # TODO: the amount should be in the model generator and be per compartments
                    'amount': self.enzyme_amounts[eid],
                    'handled_reactions': [rid]
                }
            else:
                self.enzymes[eid]['handled_reactions'].append(rid)

    def get_reaction_set_ids(self, cid, reaction_set):
        """
        :param cid: The compartment ID where the reaction set belongs
        :type cid: str
        :type reaction_set: str
        :param reaction_set:  The reaction set name of the reactions to be retrieved
        :return: A list with all the reaction ids from the specified reaction set
        within the compartment with id cid
        :rtype: list[str]
        """

        location = Location(cid, reaction_set)
        return [r for r, p in self.reactions.items() if p['location'] == location]

    def parse_reactions(self):
        """
        Parses all the reactions information and construct the reaction parameters
        for the model.

        It also generates the enzyme set obtained from reactions gene association
        """
        print '[Parser] Start parsing reactions.'
        self.enzymes = {}
        self.reactions = {}
        xml_reactions = self.model.findAll('reaction')

        bar = progressbar.ProgressBar(maxval=len(xml_reactions))
        bar.start()
        done = 0
        for reaction in xml_reactions:
            self.parse_reaction(reaction)
            bar.update(done)
            done += 1

        print '[Parser] End parsing reactions.'

    def get_compartment(self, specie_id):
        """
        :param specie_id: the specie id
        :type: str
        :return: The cid (Compartment ID) from where the specie id passed as parameter belongs
        """

        return self.model.find('species', {'id': specie_id}).get('compartment');

    def get_reaction_routing_table(self, reaction, comp_id):
        """
        :param reaction: The reaction node from the SBML xml file
        :type reaction: bs4.element.Tag
        :param comp_id: The compartment id where the reaction belongs
        :return: The routing table to know by which port should each specie be send after the
        reaction takes place and when metabolites are rejected
        :rtype: dict[str, int]
        """

        species = [s.get('species') for s in reaction.findAll('speciesReference')]
        return {s: 0 if comp_id == self.get_compartment(s)
                else 1 if self.cytoplasm_id == self.get_compartment(s)
                else 2 for s in species}

    def get_location(self, species):
        """
        :param species: A list of species ids
        :type species: list[str]
        :return: The location of the reactions conformed by the species passed as parameters
        :rtype: SBMLParser.Location
        """
        # TODO: Get a more general and flexible implementation for this in order to
        # accept different structures. An option is to load a map of str(set<compartment>) as the
        # key and the location as the value. The map can be loaded from a file allowing different
        # compartment combination interpretations.

        compartment_ids = set([self.get_compartment(s) for s in set(species)])
        location = None
        if len(compartment_ids) == 1:
            location = Location(compartment_ids.pop(), BULK)

        elif len(compartment_ids) == 2:
            if self.periplasm_id in compartment_ids:  # Periplasm inner or outer membrane

                compartment_ids.discard(self.periplasm_id)
                reaction_sets = {self.extra_cellular_id: OUTER, self.cytoplasm_id: INNER}
                location = Location(self.periplasm_id, reaction_sets[compartment_ids.pop()])
            elif self.cytoplasm_id in compartment_ids:

                compartment_ids.discard(self.cytoplasm_id)
                if self.extra_cellular_id in compartment_ids:  # Periplasm trans membrane
                    location = Location(self.periplasm_id, TRANS)
                else:  # Organelle membrane to cytoplasm
                    location = Location(compartment_ids.pop(), MEMBRANE)

        elif len(compartment_ids) == 3:  # Periplasm trans membrane
            if {self.extra_cellular_id, self.periplasm_id, self.cytoplasm_id} == compartment_ids:
                location = Location(self.periplasm_id, TRANS)

        if location is None:
            raise IllegalCompartmentCombination('Illegal compartment combination ' +
                                                str(compartment_ids))
        
        return location

    def parse_reaction(self, reaction):
        """
        Parses the reaction passed as parameter and saves in the instance attributes reactions and
        enzymes the parsed reaction parameters and the enzymes that handles the parsed reaction

        :param reaction: The reaction node from the SBML xml file
        :type: bs4.element.Tag.
        """

        rid = reaction.get('id')
        if self.is_biomass(rid) is True:
            return

        stoichiometry = self.parse_stoichiometry(reaction)
        species = stoichiometry['listOfReactants'].keys() + stoichiometry['listOfProducts'].keys()
        product_by_compartment = self.separate_by_compartment(stoichiometry['listOfProducts'])
        reactant_by_compartment = self.separate_by_compartment(stoichiometry['listOfReactants'])
        species = list(set(species))  # remove duplicates
        location = self.get_location(species)
        routing_table = self.get_reaction_routing_table(reaction, location.cid)
        parameters = {
            'rid': rid,
            'location': location,  # used by the space
            'reversible': False if reaction.get('reversible') == 'false' else True,
            'species': species,
            'product_by_compartment': product_by_compartment,
            'reactant_by_compartment': reactant_by_compartment,
            'routing_table': routing_table,
            'konSTP': self.konSTPs[rid],  # used by the space
            'konPTS': self.konPTSs[rid],  # used by the space
            'koffSTP': self.koffSTPs[rid],
            'koffPTS': self.koffPTSs[rid],
            'rate': self.rates[rid],
            'rejectTime': self.reject_times[rid]
        }

        self.parse_enzymes(reaction)
        self.reactions[rid] = parameters

    def separate_by_compartment(self, stoichiometry):
        separated_stoichiometry = {}

        for specie, amount in stoichiometry.iteritems():
            cid = self.get_compartment(specie)
            if cid not in separated_stoichiometry.keys():
                separated_stoichiometry[cid] = {specie: amount}
            else:
                separated_stoichiometry[cid][specie] = amount

        return separated_stoichiometry

    @staticmethod
    def parse_stoichiometry(reaction):
        """
        :param reaction: The reaction node from the SBML xml file
        :type: bs4.element.Tag.
        :return: The parsed stoichiometry values
        :rtype: dict[str, dict[str, float]
        """

        stoichiometry = {}

        # parsing reactants
        stoichiometries = reaction.find('listOfReactants')
        if stoichiometries is not None:
            reactants = {
                s.get('species'): 1.0
                if s.get('stoichiometry') is None
                else float(s.get('stoichiometry'))
                for s in stoichiometries.findAll('speciesReference')
            }
        else:
            reactants = {}

        # parsing products
        stoichiometries = reaction.find('listOfProducts')
        if stoichiometries is not None:
            products = {
                s.get('species'): 1.0
                if s.get('stoichiometry') is None
                else float(s.get('stoichiometry'))
                for s in stoichiometries.findAll('speciesReference')
            }
        else:
            products = {}

        stoichiometry['listOfReactants'] = reactants
        stoichiometry['listOfProducts'] = products
        return stoichiometry

    @staticmethod
    def extract_gene_association(reaction):
        """
        Extracts the GENE_ASSOCIATION: texts information from a reaction.
        If there is no associated enzymes, it returns an empty list, if there
        is some associated enzymes, it returns a list of all of them.

        :param reaction: The reaction node from the SBML xml file
        :type: bs4.element.Tag.
        :return: A list of logical gene associations.
        :rtype: list[str]
        """
        return [x.text.replace('GENE_ASSOCIATION: ', '')
                for x in reaction.notes.body.findAll('p')
                if 'GENE_ASSOCIATION' in x.text and any(char.isdigit() for char in x.text)]

    @staticmethod
    def is_biomass(rid):
        """
        :param rid: The reaction id to check.
        :type rid: str
        :return: whether the reaction id rid is a biomass rid or not
        :rtype: bool
        """
        return 'biomass' in rid.lower()
