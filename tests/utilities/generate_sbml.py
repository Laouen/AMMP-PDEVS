#!/usr/bin/python
# -*- coding: utf-8 -*-

import SBMLGenerator

def main():
reaction_amounts = {
    ('c', 'bulk'): 10,
    ('e', 'bulk'): 10,
    ('p', 'bulk'): 10,
    ('p', 'inner'): 10,
    ('p', 'outer'): 10,
    ('p', 'trans'): 10}

compartments = ['c', 'e', 'p']
species_amount = {cid: 100 for cid in compartments}

generator = SBMLGenerator.SBMLGenerator('test.xml', 'test', compartments, reaction_amounts, species_amount)
generator.create_listOfCompartments()
generator.create_listOfSpecies()
generator.create_listOfReactions()