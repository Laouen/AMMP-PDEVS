# -*- coding: utf-8 -*-

from lxml import etree
from random import randint, choice


class SBMLTestGenerator():

    def __init__(self, output, model_id, compartments, reaction_amounts, stoichimetry_elements):
        self.output_file = open(output, 'w+')
        self.compartments = compartments
        self.reaction_amounts = reaction_amounts
        self.stoichimetry_elements = stoichimetry_elements
        self.current_enzyme_id = 0
        self.species_amounts = {c: 0 for c in compartments}

        self.sbml = etree.Element('sbml')
        self.sbml.set('xmlns', 'http://www.sbml.org/sbml/level2')
        self.sbml.set('level', '2')
        self.sbml.set('version', '1')

        self.model = etree.Element('model')
        self.model.set('id', model_id)
        self.sbml.append(self.model)

        self.create_listOfCompartments()
        self.create_listOfReactions()
        self.create_listOfSpecies()

    def create_listOfCompartments(self):
        self.listOfCompartments = etree.Element('listOfCompartments')

        for cid in self.compartments:
            self.listOfCompartments.append(self.create_compartment(cid, cid + '_name'))

    def create_listOfSpecies(self):
        self.listOfSpecies = etree.Element('listOfSpecies')

        for cid, amount in self.species_amounts.items():
            for sid in (str(i) for i in range(amount)):
                self.listOfSpecies.append(self.create_specie(sid, cid))

    def create_listOfReactions(self):
        self.listOfReactions = etree.Element('listOfReactions')

        for (cid, esn), amount in self.reaction_amounts.items():
            for i in range(amount):
                rid = '_'.join(['R', str(i), cid, esn])
                self.listOfReactions.append(self.create_reaction(rid, 'false', cid, esn))

    def create_compartment(self, cid, name):
        compartment = etree.Element('compartment')
        compartment.set('id', cid)
        compartment.set('name', name)
        return compartment

    def create_specie(self, sid, compartment):
        species = etree.Element('species')
        species.set('id', sid + '_' + compartment)
        species.set('name', sid + '_name_' + compartment)
        species.set('compartment', compartment)
        return species

    def create_reaction(self, rid, reversible, cid, esn):
        reaction = etree.Element('reaction')
        reaction.set('id', rid)
        reaction.set('name', rid + '_name')
        reaction.set('reversible', reversible)

        reaction.append(self.create_note(rid))

        # the reaction has the same reactant and product
        stoichiometry = self.random_stoichiometry(cid, esn)

        listOfReactants = etree.Element('listOfReactants')
        for reactant in stoichiometry:
            listOfReactants.append(reactant)
        reaction.append(listOfReactants)

        listOfProducts = etree.Element('listOfProducts')
        for product in stoichiometry:
            listOfProducts.append(reactant)
        reaction.append(listOfProducts)

        return reaction

    def create_note(self, rid):
        note = etree.Element('notes')
        body = etree.Element('body')
        body.set('xmlns', 'http://www.w3.org/1999/xhtml')
        gene_asosiation = etree.Element('p')
        gene_asosiation.text = 'GENE_ASSOCIATION: ' + self.get_enzyme_id()
        body.append(gene_asosiation)
        note.append(body)
        return note

    def get_enzyme_id(self):
        eid = 'e' + str(self.current_enzyme_id)
        self.current_enzyme_id += 1
        return eid

    def random_stoichiometry(self, cid, esn, can_be_zero=True):

        related_compartments = [cid]
        if esn == 'outer':
            related_compartments.append('e')
        elif esn == 'inner':
            related_compartments.append('c')
        elif esn == 'trans':
            related_compartments.append('c')
            related_compartments.append('e')

        related_compartments = list(set(related_compartments))  # remove duplicated values

        for compartment in related_compartments:
            for i in range(self.stoichimetry_elements):
                s = str(self.species_amounts[compartment]) + '_' + compartment
                self.species_amounts[compartment] += 1
                specie = etree.Element('speciesReference')
                specie.set('species', s)
                specie.set('stoichiometry', '1')
                yield specie

    def save_xml(self):
        self.model.append(self.listOfCompartments)
        self.model.append(self.listOfSpecies)
        self.model.append(self.listOfReactions)

        self.output_file.write(etree.tostring(self.sbml, xml_declaration=True, encoding='UTF-8', pretty_print=True).decode('utf-8'))
        self.output_file.flush()
