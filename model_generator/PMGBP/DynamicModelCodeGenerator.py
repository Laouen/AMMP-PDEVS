#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import json
import pkg_resources

ATOMIC_MODEL_DEFINITION = pkg_resources.resource_filename(__name__, 'templates/atomic_model_definition.tpl.hpp')
DYNAMIC_ATOMIC = pkg_resources.resource_filename(__name__, 'templates/dynamic_atomic.tpl.hpp')
DYNAMIC_COUPLED = pkg_resources.resource_filename(__name__, 'templates/dynamic_coupled.tpl.hpp')
DYNAMIC_DEFINED_ATOMIC = pkg_resources.resource_filename(__name__, 'templates/dynamic_defined_atomic.tpl.hpp')
DYNAMIC_REACTION_SET = pkg_resources.resource_filename(__name__, 'templates/dynamic_reaction_set.tpl.hpp')

class DynamicModelCodeGenerator:
    """
    Writes c++ cadmium model code to the <model_name>.hpp file.
    """

    def __init__(self, model_dir='..', model_name='top', template_folder='templates', TIME='NDTime'):
        self.model_name = model_name


        # atomic ports definition template
        atomic_model_definition_tpl_path = os.curdir + os.sep + template_folder + os.sep + 'atomic_model_definition.tpl.hpp'
        self.atomic_model_def_template = open(ATOMIC_MODEL_DEFINITION, 'r').read()

        # atomic template
        atomic_tpl_path = os.curdir + os.sep + template_folder + os.sep + 'dynamic_atomic.tpl.hpp'
        self.atomic_template = open(DYNAMIC_ATOMIC, 'r').read().format(TIME=TIME)

        # coupled template
        coupled_tpl_path = os.curdir + os.sep + template_folder + os.sep + 'dynamic_coupled.tpl.hpp'
        self.coupled_template = open(DYNAMIC_COUPLED, 'r').read().format(TIME=TIME)

        # defined atomic template
        defined_atomic_tpl_path = os.curdir + os.sep + template_folder + os.sep + 'dynamic_defined_atomic.tpl.hpp'
        self.defined_atomic_template = open(DYNAMIC_DEFINED_ATOMIC, 'r').read().format(TIME=TIME)

        # reaction set definition template
        reaction_set_tpl_path = os.curdir + os.sep + template_folder + os.sep + 'dynamic_reaction_set.tpl.hpp'
        self.reaction_set_template = open(DYNAMIC_REACTION_SET, 'r').read().format(TIME=TIME)

        # port template
        self.port_template = 'struct {out_in}_{port_number}: public ' \
                             'cadmium::{out_in}_port<{message_type}>{{}};'

        # eic link template
        self.eic_template = 'cadmium::dynamic::translate::make_EIC<' \
                            '{ports_model}_ports::in_{port_number_model},' \
                            '{ports_sub_model}_ports::in_{port_number_sub_model}>' \
                            '("{id_sub_model}")'

        # eoc link template
        self.eoc_template = 'cadmium::dynamic::translate::make_EOC<' \
                            '{ports_sub_model}_ports::out_{port_number_sub_model},' \
                            '{ports_model}_ports::out_{port_number_model}>' \
                            '("{id_sub_model}")'

        # ic link template
        self.ic_template = 'cadmium::dynamic::translate::make_IC<' \
                           '{ports_model_1}_ports::out_{port_number_1},' \
                           '{ports_model_2}_ports::in_{port_number_2}>' \
                           '("{id_model_1}", "{id_model_2}")'

        self.model_file = open(model_dir + os.sep + model_name + '.hpp', 'w')
        self.model_definitions_file = open(model_dir + os.sep + model_name + '_model_definitions.hpp', 'w')

        self.tabs = ''
        self.write_includes()
        self.write_function_declaration(TIME)


    def write_reaction_set(self, cid, rsn, reaction_ids):

        reaction_ids = '{ {"' + '"}, {"'.join(['", "'.join(ids) for ids in reaction_ids]) + '"} }'

        self.write(
            self.reaction_set_template.format(
                cid=cid,
                rsn=rsn,
                reaction_ids=reaction_ids
            )
        )

        return '_'.join([cid, rsn]), 1, 3 # ioport amount of reaction sets are constants and equal to a reaction   

    def write_atomic_model(self, model_class, model_id, parameters, out_ports, in_ports, output_type, input_type):
        # ARGS mut be generated before the parameter input is modified. 
        ARGS = ', '.join(['const char*' for i in range(len(parameters))])
        ARGS = 'const char*, ' + ARGS
        parameters = ',\n\t\t'.join(map(json.dumps, parameters))
        parameters = 'xml_parameter_path.c_str(),\n\t\t' + parameters
        model_name = model_class + '_' + model_id

        self.define_atomic_model_with_ports(model_name, model_class, out_ports, in_ports, output_type, input_type)
        self.write(self.atomic_template.format(model_name=model_name, ARGS=ARGS, parameters=parameters))

        return model_name

    def write_coupled_model(self, model_id, sub_models, ports, eic, eoc, ic, bulk_solution=False):

        model_name = 'coupled_' + model_id

        ports_prefix = model_name + '_ports::'

        oiports_cpp = {'out': [], 'in': []}
        ports_cpp = []
        eic_cpp = []
        eoc_cpp = []
        ic_cpp = []

        for (port_number, message_type, out_in) in ports:
            ports_cpp.append(self.port_template.format(port_number=str(port_number),
                                                       message_type='pmgbp::types::' + message_type, # message_type namespace is not specified
                                                       out_in=out_in))
            port_name = ports_prefix + '_'.join([out_in, str(port_number)])
            oiports_cpp[out_in].append('typeid(' + port_name + ')')

        for (id_sub_model, ports_sub_model, port_number_sub_model, port_number_model) in eic:
            eic_cpp.append(self.eic_template.format(id_sub_model=id_sub_model,
                                                    ports_sub_model=ports_sub_model,
                                                    port_number_sub_model=str(port_number_sub_model),
                                                    ports_model=model_name,
                                                    port_number_model=str(port_number_model)))

        for (id_sub_model, ports_sub_model, port_number_sub_model, port_number_model) in eoc:
            eoc_cpp.append(self.eoc_template.format(id_sub_model=id_sub_model,
                                                    ports_sub_model=ports_sub_model,
                                                    port_number_sub_model=str(port_number_sub_model),
                                                    ports_model=model_name,
                                                    port_number_model=str(port_number_model)))

        for (id_model_1, ports_model_1, port_number_1, id_model_2, ports_model_2, port_number_2) in ic:
            ic_cpp.append(self.ic_template.format(id_model_1=id_model_1,
                                                  ports_model_1=ports_model_1,
                                                  port_number_1=str(port_number_1),
                                                  id_model_2=id_model_2,
                                                  ports_model_2=ports_model_2,
                                                  port_number_2=str(port_number_2)))

        self.write(self.coupled_template.format(model_name=model_name,
                                                ports='\n\t'.join(ports_cpp),
                                                oports=', '.join(oiports_cpp['out']),
                                                iports=', '.join(oiports_cpp['in']),
                                                sub_models=', '.join(sub_models),
                                                eic=', '.join(eic_cpp),
                                                eoc=', '.join(eoc_cpp),
                                                ic=', '.join(ic_cpp)))

        input_port_numbers = [port_number for (port_number, _, out_in) in ports if out_in == 'in']
        # If there is no input ports -> max = -1 -> output_ports_amount = -1 +1 = 0
        input_ports_amount = max(input_port_numbers + [-1]) + 1

        output_port_numbers = [port_number for (port_number, _, out_in) in ports if out_in == 'out']
        # If there is no output ports -> max = -1 -> output_ports_amount = -1 +1 = 0
        output_ports_amount = max(output_port_numbers + [-1]) + 1

        return model_name, input_ports_amount, output_ports_amount

    def define_atomic_model_with_ports(self, model_name, model_class, out_ports, in_ports, output_type, input_type):

        output_ports_def = '\n\t'.join([self.port_template.format(out_in='out', port_number=port_number, message_type='OUTPUT_TYPE')
                                        for port_number in range(out_ports)])

        output_port_names = ','.join(['out_' + str(port_number)
                                       for port_number in range(out_ports)])

        input_ports_def = '\n\t'.join([self.port_template.format(out_in='in', port_number=port_number, message_type='INPUT_TYPE')
                                       for port_number in range(in_ports)])

        input_port_names = ','.join(['in_' + str(port_number) for port_number in range(in_ports)])

        self.write_to_model_def(self.atomic_model_def_template.format(model_name=model_name,
                                                                      model_class=model_class,
                                                                      input_ports_definitions=input_ports_def,
                                                                      output_ports_definitions=output_ports_def,
                                                                      input_port_names=input_port_names,
                                                                      output_port_names=output_port_names,
                                                                      input_type=input_type,
                                                                      output_type=output_type))

    def write(self, code):
        code = self.tabs + code.replace('\n','\n' + self.tabs).replace('\t','    ')
        self.model_file.write(code.replace('\t','    ') + '\n')
        self.model_file.flush()

    def write_to_model_def(self, code):
        self.model_definitions_file.write(code + '\n')
        self.model_definitions_file.flush()

    def endl(self):
        self.model_file.write('\n')

    def write_includes(self):
        # model def includes
        
        structures = ['reaction', 'router', 'space']

        self.write_to_model_def('/* structure includes */\n')
        for structure in structures:
            self.write_to_model_def('#include <pmgbp/structures/' + structure + '.hpp>')
        self.write_to_model_def("")
        
        self.write('/* atomic model includes */')
        models = ['reaction', 'router', 'space']
        for model in models:
            self.write_to_model_def('#include <pmgbp/atomics/' + model + '.hpp>')
        
        self.write_to_model_def('\n#include <cadmium/modeling/ports.hpp>')

        self.write('#include <typeinfo>')
        self.write("")
        
        self.write('#include <pmgbp/model_generator/reaction_set.hpp>')

        self.write('\n#include \"' + self.model_name + '_model_definitions.hpp\"')

        self.write('#include <cadmium/modeling/dynamic_coupled.hpp>')
        self.write('#include <cadmium/modeling/dynamic_model_translator.hpp>\n\n')

    def write_function_declaration(self, TIME):
        self.write('std::shared_ptr<cadmium::dynamic::modeling::coupled<' + TIME + '>> generate_model(std::string xml_parameter_path) {')
        self.write('')
        self.tabs = '    '

    def end_model(self, top='cell'):
        self.write('return coupled_' + top + ';')
        self.tabs = ''
        self.write('}')
        self.model_file.close()
        self.model_definitions_file.close()

    # WARNING: deprecated
    def write_reaction_group(self, group_id, reaction_ids, parameters_xml):

        reaction_ids = '{ "' + '", "'.join(reaction_ids) + '" }' # vector of reaction id strings

        self.write(
            self.reaction_group_template.format(
                group_id=group_id,
                reaction_ids=reaction_ids,
                parameters_xml=parameters_xml
            )
        )

        return group_id, 1, 3 # ioport amount of reaction groups are constants 