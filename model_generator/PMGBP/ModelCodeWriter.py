#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import json
import pkg_resources

ATOMIC_MODEL_DEFINITION = pkg_resources.resource_filename(
    __name__,
    'templates/atomic_model_definition.tpl.hpp'
)
SPACE_ATOMIC_MODEL_DEFINITION = pkg_resources.resource_filename(
    __name__,
    'templates/space_atomic_model_definition.tpl.hpp'
)
ATOMIC = pkg_resources.resource_filename(
    __name__,
    'templates/atomic.tpl.hpp'
)
COUPLED = pkg_resources.resource_filename(
    __name__,
    'templates/coupled.tpl.hpp'
)
ENZYME_SET = pkg_resources.resource_filename(
    __name__,
    'templates/enzyme_set.tpl.hpp'
)


class ModelCodeWriter:
    """
    Writes c++ cadmium model code to the <model_name>.hpp file.
    """

    def __init__(self, model_dir='..', model_name='top', time='NDTime'):
        self.model_name = model_name

        # space atomic model definition template
        self.space_atomic_model_def_template = open(SPACE_ATOMIC_MODEL_DEFINITION, 'r').read()

        # atomic template
        self.atomic_template = open(ATOMIC, 'r').read().format(TIME=time)

        # coupled template
        self.coupled_template = open(COUPLED, 'r').read().format(TIME=time)

        # enzyme set definition template
        self.enzyme_set_template = open(ENZYME_SET, 'r').read().format(TIME=time)

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
        self.write_function_declaration(time)

    def write_space_atomic_model(
            self,
            model_id,
            parameters,
            out_ports,
            in_ports,
            product_type,
            reactant_type,
            information_type
    ):
        # ARGS mut be generated before the parameter input is modified. 
        args = ', '.join(['const char*' for _ in range(len(parameters))])
        args = 'const char*, ' + args
        parameters = ',\n\t\t'.join(map(json.dumps, parameters))
        parameters = 'xml_parameter_path.c_str(),\n\t\t' + parameters
        model_name = 'space_' + model_id

        self.define_space_atomic_model_with_ports(
            model_name,
            out_ports,
            in_ports,
            product_type,
            reactant_type,
            information_type
        )
        self.write(self.atomic_template.format(model_name=model_name, ARGS=args, parameters=parameters))

        return model_name

    def define_space_atomic_model_with_ports(
            self,
            model_name,
            out_ports,
            in_ports,
            product_type,
            reactant_type,
            information_type
    ):

        # Create output ports
        output_ports_def = [
            self.port_template.format(
                out_in='out',
                port_number=port_number,
                message_type='pmgbp::types::' + reactant_type
            ) for port_number in range(out_ports)
        ]

        output_ports_def = '\n\t'.join(output_ports_def)

        output_port_names = ['out_' + str(port_number) for port_number in range(out_ports)]
        output_port_names = ','.join(output_port_names)

        # Create input ports
        # Input ports are duplacated to recive product and information messages
        input_port_defs = [
            self.port_template.format(
                out_in='in',
                port_number=str(port_number) + '_product',
                message_type='pmgbp::types::' + product_type) for port_number in range(in_ports)
        ]
        input_port_defs += [
            self.port_template.format(
                out_in='in',
                port_number=str(port_number) + '_information',
                message_type='pmgbp::types::' + information_type) for port_number in range(in_ports)
        ]

        input_ports_def = '\n\t'.join(input_port_defs)

        input_port_names = ['in_' + str(port_number) + '_product' for port_number in range(in_ports)]
        input_port_names += ['in_' + str(port_number) + '_information' for port_number in range(in_ports)]
        input_port_names = ','.join(input_port_names)

        self.write_to_model_def(
            self.space_atomic_model_def_template.format(
                model_name=model_name,
                input_ports_definitions=input_ports_def,
                output_ports_definitions=output_ports_def,
                input_port_names=input_port_names,
                output_port_names=output_port_names,
                product_type=product_type,
                reactant_type=reactant_type,
                information_type=information_type
            )
        )

    def write_enzyme_set(self, cid, esn, enzyme_ids):
        enzyme_ids = '{ {"' + '"}, {"'.join(['", "'.join(ids) for ids in enzyme_ids]) + '"} }'

        self.write(
            self.enzyme_set_template.format(
                cid=cid,
                esn=esn,
                enzyme_ids=enzyme_ids
            )
        )

        return '_'.join([cid, esn]), 1, 3  # ioport amount of reaction sets are constants and equal to a reaction

    def write_coupled_model(self, model_id, sub_models, ports, eic, eoc, ic):

        model_name = 'coupled_' + model_id

        ports_prefix = model_name + '_ports::'

        oiports_cpp = {'out': [], 'in': []}
        ports_cpp = []
        eic_cpp = []
        eoc_cpp = []
        ic_cpp = []

        for (port_number, message_type, out_in) in ports:
            ports_cpp.append(
                self.port_template.format(
                    port_number=str(port_number),
                    message_type='pmgbp::types::' + message_type,  # message_type namespace is not specified
                    out_in=out_in
                )
            )
            port_name = ports_prefix + '_'.join([out_in, str(port_number)])
            oiports_cpp[out_in].append('typeid(' + port_name + ')')

        for (id_sub_model, ports_sub_model, port_number_sub_model, port_number_model) in eic:
            eic_cpp.append(
                self.eic_template.format(
                    id_sub_model=id_sub_model,
                    ports_sub_model=ports_sub_model,
                    port_number_sub_model=str(port_number_sub_model),
                    ports_model=model_name,
                    port_number_model=str(port_number_model)
                )
            )

        for (id_sub_model, ports_sub_model, port_number_sub_model, port_number_model) in eoc:
            eoc_cpp.append(
                self.eoc_template.format(
                    id_sub_model=id_sub_model,
                    ports_sub_model=ports_sub_model,
                    port_number_sub_model=str(port_number_sub_model),
                    ports_model=model_name,
                    port_number_model=str(port_number_model)
                )
            )

        for (id_model_1, ports_model_1, port_number_1, id_model_2, ports_model_2, port_number_2) in ic:
            ic_cpp.append(
                self.ic_template.format(
                    id_model_1=id_model_1,
                    ports_model_1=ports_model_1,
                    port_number_1=str(port_number_1),
                    id_model_2=id_model_2,
                    ports_model_2=ports_model_2,
                    port_number_2=str(port_number_2)
                )
            )

        ports_str = '' if len(ports_cpp) == 0 else '\n\t' + '\n\t'.join(ports_cpp) + '\n'
        oports_str = '' if len(oiports_cpp['out']) == 0 else '\n\t' + ',\n\t'.join(oiports_cpp['out']) + '\n'
        iports_str = '' if len(oiports_cpp['in']) == 0 else '\n\t' + ',\n\t'.join(oiports_cpp['in']) + '\n'
        sub_models_str = '' if len(sub_models) == 0 else '\n\t' + ',\n\t'.join(sub_models) + '\n'
        eic_str = '' if len(eic_cpp) == 0 else '\n\t' + ',\n\t'.join(eic_cpp) + '\n'
        eoc_str = '' if len(eoc_cpp) == 0 else '\n\t' + ',\n\t'.join(eoc_cpp) + '\n'
        ic_str = '' if len(ic_cpp) == 0 else '\n\t' + ',\n\t'.join(ic_cpp) + '\n'

        self.write(
            self.coupled_template.format(
                model_name=model_name,
                ports=ports_str,
                oports=oports_str,
                iports=iports_str,
                sub_models=sub_models_str,
                eic=eic_str,
                eoc=eoc_str,
                ic=ic_str
            )
        )

        input_port_numbers = [
            int(str(port_number).split('_')[0])
            for (port_number, _, out_in) in ports if out_in == 'in'
        ]
        # If there is no input ports -> max = -1 -> output_ports_amount = -1 +1 = 0
        input_ports_amount = max(input_port_numbers + [-1]) + 1

        output_port_numbers = [
            int(str(port_number).split('_')[0])
            for (port_number, _, out_in) in ports if out_in == 'out'
        ]
        # If there is no output ports -> max = -1 -> output_ports_amount = -1 +1 = 0
        output_ports_amount = max(output_port_numbers + [-1]) + 1

        return model_name, input_ports_amount, output_ports_amount

    def write(self, code):
        code = self.tabs + code.replace('\n', '\n' + self.tabs).replace('\t', '    ')
        self.model_file.write(code.replace('\t', '    ') + '\n')
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
        models = ['enzyme', 'router', 'space']
        for model in models:
            self.write_to_model_def('#include <pmgbp/atomics/' + model + '.hpp>')
        
        self.write_to_model_def('\n#include <cadmium/modeling/ports.hpp>')

        self.write('#include <typeinfo>')
        self.write("")
        
        self.write('#include <pmgbp/model_generator/enzyme_set.hpp>')

        self.write('\n#include \"' + self.model_name + '_model_definitions.hpp\"')

        self.write('#include <cadmium/modeling/dynamic_coupled.hpp>')
        self.write('#include <cadmium/modeling/dynamic_model_translator.hpp>\n\n')

    def write_function_declaration(self, time):
        self.write(
            ('std::shared_ptr<cadmium::dynamic::modeling::coupled<'
             + time
             + '>> generate_model(std::string xml_parameter_path) {')
        )
        self.write('')
        self.tabs = '    '

    def end_model(self, top='cell'):
        self.write('return coupled_' + top + ';')
        self.tabs = ''
        self.write('}')
        self.model_file.close()
        self.model_definitions_file.close()
