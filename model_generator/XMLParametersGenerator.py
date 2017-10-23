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

        router = etree.Element(model_id)
        routing_table = self.generate_table(routing_table,
                                            'routingTable',
                                            'entry',
                                            ['metaboliteId'],
                                            'port')
        router.append(routing_table)
        self.routers.append(router)

    def add_reaction(self, model_id, parameters):

        reaction = etree.Element(model_id)

        parameter_keys = ['rate', 'rejectTime', 'koffSTP', 'koffPTS']
        for key in parameter_keys:
            parameter = etree.Element(key)
            parameter.text = str(parameters[key])
            reaction.append(parameter)

        routing_table = self.generate_table(parameter['routing_table'],
                                            'routingTable',
                                            'entry',
                                            ['metaboliteId'],
                                            'port')
        reaction.append(routing_table)

        stoichiometry_by_compartments = etree.Element('stoichiometryByCompartments')

        product_compartments = parameters.product_by_compartment.keys()
        reactant_compartments = parameters.reactant_by_compartment.keys()

        compartments = set(product_compartments + reactant_compartments)

        for cid in compartments:
            compartment = etree.Element('compartment')
            xml_cid = etree.Element('id')
            xml_cid.text = cid
            compartment.append(xml_cid)

            if cid in reactant_compartments:
                stoichiometry = parameters.product_by_compartment[cid]
                stoichiometry = self.generate_table(stoichiometry,
                                                    'substrate',
                                                    'specie',
                                                    ['id'],
                                                    'amount')
                compartment.append(stoichiometry)

            if cid in product_compartments:
                stoichiometry = parameters.product_by_compartment[cid]
                stoichiometry = self.generate_table(stoichiometry,
                                                    'product',
                                                    'specie',
                                                    ['id'],
                                                    'amount')
                compartment.append(stoichiometry)

            stoichiometry_by_compartments.append(compartment)

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
