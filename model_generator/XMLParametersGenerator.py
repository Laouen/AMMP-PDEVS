#!/usr/bin/python
# -*- coding: utf-8 -*-

from lxml import etree


class XMLParametersGenerator:

    def __init__(self):
        self.parameters = etree.Element('parameters')

        self.spaces = etree.Element('spaces')
        self.routers = etree.Element('routers')
        self.reactions = etree.Element('reactions')

        self.parameters.append(self.spaces)
        self.parameters.append(self.routers)
        self.parameters.append(self.reactions)

    def add_router(self, model_id, routing_table):

        xml_router = etree.Element(model_id)
        xml_routing_table = self.generate_table(routing_table,
                                                'routingTable',
                                                'entry',
                                                ['metaboliteId'],
                                                'port')
        xml_router.append(xml_routing_table)
        self.routers.append(xml_router)

    def add_reaction(self, model_id, parameters):

        xml_reaction = etree.Element(model_id)

        parameter_keys = ['rate', 'rejectTime', 'koffSTP', 'koffPTS']
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

        product_compartments = parameters.product_by_compartment.keys()
        reactant_compartments = parameters.reactant_by_compartment.keys()

        compartments = set(product_compartments + reactant_compartments)

        for cid in compartments:
            xml_compartment = etree.Element('compartment')

            xml_cid = etree.Element('id')
            xml_cid.text = cid
            xml_compartment.append(xml_cid)

            if cid in reactant_compartments:
                stoichiometry = parameters.product_by_compartment[cid]
                xml_stoichiometry = self.generate_table(stoichiometry,
                                                        'substrate',
                                                        'specie',
                                                        ['id'],
                                                        'amount')
                xml_compartment.append(xml_stoichiometry)

            if cid in product_compartments:
                stoichiometry = parameters.product_by_compartment[cid]
                xml_stoichiometry = self.generate_table(stoichiometry,
                                                        'product',
                                                        'specie',
                                                        ['id'],
                                                        'amount')
                xml_compartment.append(xml_stoichiometry)

            xml_stoichiometry_by_compartments.append(xml_compartment)
        xml_reaction.append(xml_stoichiometry_by_compartments)

    def add_space(self, model_id, parameters, routing_table):
        xml_space = etree.Element(model_id)

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
                                                ['cid', 'rsn'],
                                                'port')
        xml_space.append(xml_routing_table)

        cid = parameters['cid']
        xml_enzymes = etree.Element('enzymes')
        for eid, enzyme_parameters in parameters['enzymes']:
            xml_enzyme = etree.Element('enzyme')

            parameter_keys = ['id', 'amount']
            for key in parameter_keys:
                xml_parameter = etree.Element(key)
                xml_parameter.text = str(enzyme_parameters[key])
                xml_enzyme.append(xml_parameter)

            for rid in enzyme_parameters['handled_reactions']:
                reaction_parameters = parameters['reaction_parameters'][rid]

                xml_reaction = etree.Element('reaction')

                parameter_keys = ['rid', 'rejectTime', 'konSTP', 'konPTS',
                                  'koffSTP', 'koffPTS', 'reversible']
                for key in parameter_keys:
                    parameter = etree.Element(key)
                    parameter.text = str(reaction_parameters[key])
                    xml_reaction.append(parameter)

                xml_address = etree.Element('address')
                xml_address.set('cid', reaction_parameters['location'].cid)
                xml_address.set('rsn', reaction_parameters['location'].rsn)
                xml_reaction.append(xml_address)

                stoichiometry = etree.Element('stoichiometry')

                product = reaction_parameters['product_by_compartment'][cid]
                xml_product = self.generate_table(product,
                                                  'product',
                                                  'specie',
                                                  ['id'],
                                                  'amount')
                stoichiometry.append(xml_product)

                substrate = reaction_parameters['reactant_by_compartment'][cid]
                xml_substrate = self.generate_table(substrate,
                                                    'substrate',
                                                    'specie',
                                                    ['id'],
                                                    'amount')
                stoichiometry.append(xml_substrate)
            xml_enzymes.append(xml_enzyme)
        xml_space.append(xml_enzymes)

    def print_parameters(self):
        print etree.tostring(self.parameters, encoding='UTF-8', pretty_print=True)

    @staticmethod
    def generate_table(table, table_name, entry_name, attributes, value_attribute):
        xml_table = etree.Element(table_name)

        for key_values, port_number in table.iteritems():
            entry = etree.Element(entry_name)

            if type(key_values) is not tuple:
                key_values = tuple([key_values])

            for k, v in zip(attributes, key_values):
                entry.set(k, str(v))

            entry.set(value_attribute, str(port_number))
            xml_table.append(entry)

        return xml_table
