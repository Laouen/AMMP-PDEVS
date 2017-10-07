# -*- coding: utf-8 -*-
import os

class ModelCodeGenerator:
    """
    Writes c++ cadmium model code to the <model_name>.hpp file.
    """

    def __init__(self, model_name='top', template_folder='templates'):
        self.model_name = model_name

        # atomic template
        atomic_tpl_path = os.curdir + os.sep + template_folder + os.sep + 'atomic.tpl.hpp'
        self.atomic_template = open(atomic_tpl_path, 'r').read()

        # atomic ports template
        ports_tpl_path = os.curdir + os.sep + template_folder + os.sep + 'ports.tpl.hpp'
        self.ports_def_template = open(ports_tpl_path, 'r').read()

        # coupled template
        coupled_tpl_path = os.curdir + os.sep + template_folder + os.sep + 'coupled.tpl.hpp'
        self.coupled_template = open(coupled_tpl_path, 'r').read()

        # port template
        self.atomic_port_template = 'struct {out_in}_{port_number} : public ' \
                                    'cadmium::{out_in}_port<{message_type}> {{}};'

        self.coupled_port_template = 'struct {model_name}_{out_in}_{port_number}: public ' \
                                     'cadmium::{out_in}_port<pmgbp::types::{message_type}>{{}};'

        self.eic_template = 'cadmium::modeling::EIC<' \
                            '{model_name}_in_{model_port_number},' \
                            '{sub_model_name},' \
                            '{sub_model_name}_ports::in_{sub_model_port_number}>'

        self.eoc_template = 'cadmium::modeling::EOC<' \
                            '{sub_model_name},' \
                            '{sub_model_name}_ports::out_{sub_model_port_number},' \
                            '{model_name}_out_{model_port_number}>'

        self.ic_template =  'cadmium::modeling::IC<' \
                            '{sub_model_1},{sub_model_1}_ports::out_{port_number_1},' \
                            '{sub_model_2},{sub_model_2}_ports::in_{port_number_2}>'

        self.model_file = open(model_name + '.hpp', 'wb')
        self.port_file = open(model_name + '_ports.hpp', 'wb')

    def write_atomic_model(self, model_class, model_id, parameters, out_ports, in_ports, output_type, input_type):

        model_name = model_class + "_" + model_id

        self.write_atomic_model_ports(model_name, out_ports, in_ports)

        self.write(self.atomic_template.format(model_name=model_name,
                                               model_class=model_class,
                                               parameters=', '.join([model_id] + parameters),
                                               output_type=output_type,
                                               input_type=input_type))

    def write_coupled_model(self, model_id, submodels, ports, eic, eoc, ic):

        model_name = "coupled_" + model_id

        oiports_cpp = {'out': [], 'in': []}
        ports_cpp = []
        eic_cpp = []
        eoc_cpp = []
        ic_cpp = []

        for (port_number, message_type, out_in) in ports:
            ports_cpp.append(self.coupled_port_template.format(model_name=model_name,
                                                               port_number=str(port_number),
                                                               message_type=message_type,
                                                               out_in=out_in))
            oiports_cpp[out_in].append('_'.join([model_name, out_in, str(port_number)]))

        for (model_port_number, sub_model_name, sub_model_port_number) in eic:
            eic_cpp.append(self.eic_template.format(model_name=model_name,
                                                    model_port_number=str(model_port_number),
                                                    sub_model_name=sub_model_name,
                                                    sub_model_port_number=str(sub_model_port_number)))

        for (sub_model_name, sub_model_port_number, model_port_number) in eoc:
            eoc_cpp.append(self.eoc_template.format(model_name=model_name,
                                                    sub_model_name=sub_model_name,
                                                    sub_model_port_number=str(sub_model_port_number),
                                                    model_port_number=str(model_port_number)))

        for (sub_model_1, port_number_1, sub_model_2, port_number_2) in ic:
            ic_cpp.append(self.ic_template.format(sub_model_1=sub_model_1,
                                                  port_number_1=str(port_number_1),
                                                  sub_model_2=sub_model_2,
                                                  port_number_2=str(port_number_2)))

        self.write(self.coupled_template.format(model_name=model_name,
                                                ports='\n'.join(ports_cpp),
                                                oports=', '.join(oiports_cpp['out']),
                                                iports=', '.join(oiports_cpp['in']),
                                                submodels=', '.join(submodels),
                                                eic=', '.join(eic_cpp),
                                                eoc=', '.join(eoc_cpp),
                                                ic=', '.join(ic_cpp)))

    def write_atomic_model_ports(self, model_name, out_ports, in_ports):

        output_ports_def = '\n\t'.join([self.atomic_port_template.format(out_in='out',
                                                                         port_number=port_number,
                                                                         message_type='OUTPUT_TYPE')
                                      for port_number in range(out_ports)])

        output_ports_names = ','.join(['out_' + str(port_number) for port_number in range(out_ports)])

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
