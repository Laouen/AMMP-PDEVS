# -*- coding: utf-8 -*-
import os
import json


class DynamicModelCodeGenerator:
    """
    Writes c++ cadmium model code to the <model_name>.hpp file.
    """

    def __init__(self, model_dir='..', model_name='top', template_folder='templates', TIME='NDTime'):
        self.model_name = model_name

        # defined atomic template
        defined_atomic_tpl_path = os.curdir + os.sep + template_folder + os.sep + 'dynamic_defined_atomic.tpl.hpp'
        self.defined_atomic_template = open(defined_atomic_tpl_path, 'r').read().format(TIME=TIME)

        # atomic template
        atomic_tpl_path = os.curdir + os.sep + template_folder + os.sep + 'dynamic_atomic.tpl.hpp'
        self.atomic_template = open(atomic_tpl_path, 'r').read().format(TIME=TIME)

        # atomic ports definition template
        ports_tpl_path = os.curdir + os.sep + template_folder + os.sep + 'ports.tpl.hpp'
        self.ports_def_template = open(ports_tpl_path, 'r').read()

        # coupled template
        coupled_tpl_path = os.curdir + os.sep + template_folder + os.sep + 'dynamic_coupled.tpl.hpp'
        self.coupled_template = open(coupled_tpl_path, 'r').read().format(TIME=TIME)

        # reaction group definition template
        self.reaction_group_template = 'std::shared_ptr<cadmium::dynamic::modeling::coupled<'+ TIME +'>> reaction_group_{{group_id}} =' \
                                       'make_reaction_group("{group_id}", {reaction_ids}, "{parameter_xml}");'

        # todo: unify atomic port template and coupled port template, the only differ because the atomic one 
        # defines the message type all in the template while the coupled one hardoces the pmgbp::types:: namespace
        # atomic model port template
        self.atomic_port_template = 'struct {out_in}_{port_number}: public ' \
                                    'cadmium::{out_in}_port<{message_type}>{{}};'

        # coupled model port template
        self.port_template = 'struct {out_in}_{port_number}: public ' \
                             'cadmium::{out_in}_port<pmgbp::types::{message_type}>{{}};'

        # eic link template
        self.eic_template = 'cadmium::dynamic::translate::make_EIC<' \
                            '{model_name}_ports::in_{model_port_number},' \
                            '{sub_model_name}_ports::in_{sub_model_port_number}>' \
                            '("{sub_model_id}")'

        # eoc link template
        self.eoc_template = 'cadmium::dynamic::translate::make_EOC<' \
                            '{sub_model_name}_ports::out_{sub_model_port_number},' \
                            '{model_name}_ports::out_{model_port_number}>' \
                            '("{sub_model_id}")'

        # ic link template
        self.ic_template = 'cadmium::dynamic::translate::make_IC<' \
                           '{sub_model_1}_ports::out_{port_number_1},' \
                           '{sub_model_2}_ports::in_{port_number_2}>' \
                           '("{sub_model_1_id}", "{sub_model_2_id}")'

        self.model_file = open(model_dir + os.sep + model_name + '.hpp', 'wb')
        self.port_file = open(model_dir + os.sep + model_name + '_ports.hpp', 'wb')

        self.write_includes()

    # thi method does the same as write_atomic_model but instead of declaring a new model variable
    # it pushs the new model to a std::vector of models pass as parameter
    # WARNING: this method is only made to construct reaction atomic models
    def write_reaction_group_template(self,
                                       model_ids,
                                       group_id,
                                       parameters_xml):

        self.write(self.reaction_group_template.format(group_id=group_id, parameters_xml=parameters_xml))



    def write_atomic_model(self,
                           model_class,
                           model_id,
                           parameters,
                           out_ports,
                           in_ports,
                           output_type,
                           input_type,
                           is_defined=False):

        # ARGS mut be generated before the parameter input is modified. 
        ARGS = ', '.join(['const char*' for i in range(len(parameters))])
        model_name = model_class + '_' + model_id
        parameters = ',\n\t\t'.join(map(json.dumps, parameters))

        if not is_defined:
            self.write_atomic_model_ports(model_name, out_ports, in_ports)
            self.write(self.atomic_template.format(model_name=model_name,
                                                   model_class=model_class,
                                                   output_type=output_type,
                                                   input_type=input_type,
                                                   ARGS=ARGS,
                                                   parameters=parameters))
        else:
            self.write(self.defined_atomic_template.format(model_name=model_name,
                                                           model_class=model_class,
                                                           ARGS=ARGS,
                                                           parameters=parameters))

        return model_name

    def write_coupled_model(self, model_id, sub_models, ports, eic, eoc, ic):

        model_name = 'coupled_' + model_id

        port_prefix = model_name + '_ports::'

        oiports_cpp = {'out': [], 'in': []}
        ports_cpp = []
        eic_cpp = []
        eoc_cpp = []
        ic_cpp = []

        for (port_number, message_type, out_in) in ports:
            ports_cpp.append(self.port_template.format(port_number=str(port_number),
                                                       message_type=message_type,
                                                       out_in=out_in))
            port_name = port_prefix + '_'.join([out_in, str(port_number)])
            oiports_cpp[out_in].append('typeid(' + port_name + ')')

        for (model_port_number, sub_model_name, sub_model_port_number) in eic:
            sub_model_id = sub_model_name
            if 'router' in sub_model_name:
                sub_model_name = 'pmgbp::models::router'
            eic_cpp.append(self.eic_template.format(model_name=model_name,
                                                    model_port_number=str(model_port_number),
                                                    sub_model_name=sub_model_name,
                                                    sub_model_id=sub_model_id,
                                                    sub_model_port_number=str(sub_model_port_number)))

        for (sub_model_name, sub_model_port_number, model_port_number) in eoc:
            sub_model_id = sub_model_name
            if 'reaction' in sub_model_name:
                sub_model_name = 'pmgbp::models::reaction'

            eoc_cpp.append(self.eoc_template.format(model_name=model_name,
                                                    sub_model_name=sub_model_name,
                                                    sub_model_port_number=str(sub_model_port_number),
                                                    model_port_number=str(model_port_number),
                                                    sub_model_id=sub_model_id))

        for (sub_model_1, port_number_1, sub_model_2, port_number_2) in ic:
            sub_model_1_id = sub_model_1
            sub_model_2_id = sub_model_2
            if 'reaction' in sub_model_2:
                sub_model_2 = 'pmgbp::models::reaction'
            
            if 'router' in sub_model_1:
                sub_model_1 = 'pmgbp::models::router'


            ic_cpp.append(self.ic_template.format(sub_model_1=sub_model_1,
                                                  sub_model_1_id=sub_model_1_id,
                                                  port_number_1=str(port_number_1),
                                                  sub_model_2=sub_model_2,
                                                  sub_model_2_id=sub_model_2_id,
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

    def write_atomic_model_ports(self, model_name, out_ports, in_ports):

        output_ports_def = '\n\t'.join([self.atomic_port_template.format(out_in='out',
                                                                         port_number=port_number,
                                                                         message_type='OUTPUT_TYPE')
                                        for port_number in range(out_ports)])

        output_ports_names = ','.join(['out_' + str(port_number)
                                       for port_number in range(out_ports)])

        input_ports_def = '\n\t'.join([self.atomic_port_template.format(out_in='in',
                                                                        port_number=port_number,
                                                                        message_type='INPUT_TYPE')
                                       for port_number in range(in_ports)])

        input_ports_names = ','.join(['in_' + str(port_number) for port_number in range(in_ports)])

        self.write_ports(self.ports_def_template.format(model_name=model_name,
                                                        output_ports_definitions=output_ports_def,
                                                        input_ports_definitions=input_ports_def,
                                                        input_port_names=input_ports_names,
                                                        output_port_names=output_ports_names))

    def write(self, code):
        self.model_file.write(code + '\n')
        self.model_file.flush()

    def write_ports(self, code):
        self.port_file.write(code + '\n')
        self.port_file.flush()

    def endl(self):
        self.model_file.write('\n')

    def write_includes(self):
        # model def includes
        
        self.write('#include <typeinfo>')

        structures = ['reaction', 'router', 'space']
        models = ['reaction', 'router', 'space']

        self.write_ports('/* structure includes */\n')
        for structure in structures:
            self.write_ports('#include <pmgbp/structures/' + structure + '.hpp>')
        self.write_ports("")

        self.write('\n/* atomic model includes */')
        for model in models:
            self.write('#include <pmgbp/atomics/' + model + '.hpp>')

        self.write('\n#include \"' + self.model_name + '_ports.hpp\"')

        self.write('\n#include <cadmium/modeling/ports.hpp>')
        self.write('#include <cadmium/modeling/dynamic_coupled.hpp>')
        self.write('#include <cadmium/modeling/dynamic_model_translator.hpp>\n\n')
