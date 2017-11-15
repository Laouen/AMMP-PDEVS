#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import gflags
import pandas as pd
import json


def main(FLAGS):

    current_time = ""
    data = {'time': [], 'specie': [], 'amount': []}

    with open(FLAGS.input_file, 'r') as file:

        for line in file.readlines():
            if line.startswith('State for model space'):
                get_data_from_state(data, current_time, line)
            elif not line.startswith('[') and not line.startswith('State for model'):
                current_time = line.replace('\n', '').replace('\"', '')

    df = pd.DataFrame(data)
    df.to_csv(FLAGS.output_file, index=False, encoding='utf-8', sep=';')


def get_data_from_state(data, current_time, line):
    json_data = line[line.index('is') + 3:]
    json_data = json.loads(json_data)

    for specie, amount in json_data['metabolites'].iteritems():
        data['specie'].append(specie)
        data['amount'].append(amount)
        data['time'].append(current_time)
    return 0


if __name__ == '__main__':

    gflags.DEFINE_string('output_file', 'output.csv', 'The output csv file', short_name='o')
    gflags.DEFINE_string('input_file', 'output.log', 'The model log file to extract the data', short_name='f')

    FLAGS = gflags.FLAGS

    try:
        argv = FLAGS(sys.argv)
    except gflags.FlagsError as e:
        print '%s \nUsage: %s \n%s' % (e, sys.argv[0], FLAGS)
        sys.exit(1)

    main(FLAGS)
