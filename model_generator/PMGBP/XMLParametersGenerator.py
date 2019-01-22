#!/usr/bin/python3
# -*- coding: utf-8 -*-

from lxml import etree
import os


class XMLParametersGenerator:

    def __init__(self, model_dir='..', xml_file='parameters'):

        self.xml_file = open(model_dir + os.sep + xml_file + '.xml', 'wb')
        self.xml_file_path = self.xml_file.name

        self.parameters = etree.Element('parameters')

        self.spaces = etree.Element('spaces')
        self.routers = etree.Element('routers')
        self.enzymes = etree.Element('enzymes')
        self.reactions = etree.Element('reactions')

        self.parameters.append(self.spaces)
        self.parameters.append(self.routers)
        self.parameters.append(self.enzymes)
        self.parameters.append(self.reactions)

    def add_router(self, model_id, routing_table):

        xml_router = etree.Element(model_id)
        xml_routing_table = self.generate_table(routing_table,
                                                'routingTable',
                                                'entry',
                                                ['enzymeID'],
                                                'port')
        xml_router.append(xml_routing_table)
        self.routers.append(xml_router)

    def add_reaction(self, rid, parameters):

        xml_reaction = etree.Element(rid)

        parameter_keys = ['rate', 'rejectRate', 'koffSTP', 'koffPTS', 'konSTP', 'konPTS', 'reversible']
        for key in parameter_keys:
            xml_parameter = etree.Element(key)
            xml_parameter.text = str(parameters[key])
            xml_reaction.append(xml_parameter)

        xml_routing_table = self.generate_table(parameters['routing_table'],
                                                'routingTable',
                                                'entry',
                                                ['metaboliteId'],
                                                'port')
        xml_reaction.append(xml_routing_table)

        xml_stoichiometry_by_compartments = etree.Element('stoichiometryByCompartments')

        product_compartments = list(parameters['product_by_compartment'].keys())
        reactant_compartments = list(parameters['reactant_by_compartment'].keys())

        compartments = set(product_compartments + reactant_compartments)

        for cid in compartments:
            xml_compartment = etree.Element('compartmentStoichiometry')
            xml_compartment.set('cid', cid)

            if cid in reactant_compartments:
                stoichiometry = parameters['reactant_by_compartment'][cid]
                xml_stoichiometry = self.generate_table(stoichiometry,
                                                        'substrate',
                                                        'specie',
                                                        ['id'],
                                                        'amount')
                xml_compartment.append(xml_stoichiometry)

            if cid in product_compartments:
                stoichiometry = parameters['product_by_compartment'][cid]
                xml_stoichiometry = self.generate_table(stoichiometry,
                                                        'product',
                                                        'specie',
                                                        ['id'],
                                                        'amount')
                xml_compartment.append(xml_stoichiometry)

            xml_stoichiometry_by_compartments.append(xml_compartment)
        xml_reaction.append(xml_stoichiometry_by_compartments)
        self.reactions.append(xml_reaction)

    def add_enzyme(self, model_id, parameters):

        xml_enzyme = etree.Element(model_id)
        xml_enzyme_reactions = etree.Element("reactions")

        for rid in parameters['handled_reactions']:
            xml_reaction = etree.Element("reaction")
            xml_reaction.set("id", rid)
            xml_enzyme_reactions.append(xml_reaction)

        xml_enzyme.append(xml_enzyme_reactions)

        self.enzymes.append(xml_enzyme)

    def add_space(self, model_id, parameters, routing_table):
        xml_space = etree.Element(model_id)

        xml_volume = etree.Element('volume')
        xml_volume.text = str(parameters['volume'])

        xml_space.append(xml_volume)

        xml_interval_time = etree.Element('intervalTime')
        xml_interval_time.text = parameters['interval_time']

        xml_space.append(xml_interval_time)

        xml_metabolites = self.generate_table(parameters['metabolites'],
                                              'metabolites',
                                              'metabolite',
                                              ['id'],
                                              'amount')
        xml_space.append(xml_metabolites)

        xml_routing_table = self.generate_table(routing_table,
                                                'routingTable',
                                                'entry',
                                                ['cid', 'esn'],
                                                'port')
        xml_space.append(xml_routing_table)

        cid = parameters['cid']
        xml_enzymes_parameters = etree.Element('enzymes')

        for eid, enzyme_parameters in parameters['enzymes'].items():
            
            # separate related enzymes by location
            space_related_reactions = list(parameters['reaction_parameters'].keys())
            all_locations = set([parameters['reaction_parameters'][rid]['location'] 
                                 for rid in enzyme_parameters['handled_reactions'] 
                                 if rid in space_related_reactions])
    
            for location in all_locations:
                xml_enzyme_parameters = etree.Element(enzyme_parameters['id'])
                xml_enzyme_parameters.set('amount', str(enzyme_parameters['amount']))

                # Set enzyme address
                xml_address = etree.Element('address')
                xml_address.set('cid', location.cid)
                xml_address.set('esn', location.esn)
                xml_enzyme_parameters.append(xml_address)

                # Set enzyme handled reactions in the location
                xml_handled_reactions = etree.Element('reactions')
                handled_reactions = [rid for rid in enzyme_parameters['handled_reactions']
                                     if rid in space_related_reactions
                                     and parameters['reaction_parameters'][rid]['location'] == location]

                # Set reactions' parameters 
                for rid in handled_reactions:
                    reaction_parameters = parameters['reaction_parameters'][rid]
                    xml_reaction = etree.Element('reaction')
                    xml_reaction.set('id', reaction_parameters['rid'])

                    xml_handled_reactions.append(xml_reaction)
            
                xml_enzyme_parameters.append(xml_handled_reactions)
                xml_enzymes_parameters.append(xml_enzyme_parameters)

            xml_space.append(xml_enzymes_parameters)

        self.spaces.append(xml_space)

    def print_parameters(self):
        print(etree.tostring(self.parameters, encoding='UTF-8', pretty_print=True))

    def save_xml(self):
        self.xml_file.write(etree.tostring(self.parameters, encoding='UTF-8', pretty_print=True))
        self.xml_file.flush()

    @staticmethod
    def generate_table(table, table_name, entry_name, attributes, value_attribute):
        xml_table = etree.Element(table_name)

        for key_values, port_number in table.items():
            entry = etree.Element(entry_name)

            if type(key_values) is not tuple:
                key_values = tuple([key_values])

            for k, v in zip(attributes, key_values):
                entry.set(k, str(v))

            entry.set(value_attribute, str(port_number))
            xml_table.append(entry)

        return xml_table
