#!/usr/bin/env python3

import sys
import argparse
import os
import json

parser = argparse.ArgumentParser(description='Quickly patch KF PT ID update for PNOC WES')
parser.add_argument('-d', '--data', action='store', dest='data', help='Sample file with patients IDS and external IDs')
parser.add_argument('-s', '--sample', action='store', dest='sample', help='Sample file with patients IDS and external IDs')
parser.add_argument('-j', '--json', action='store', dest='json', help='Input json file')

args = parser.parse_args()
with open(args.json) as f:
    pt_data = json.load(f)
pt_id_dict = {}
pt_sheet = open('data_clinical_patient.txt', 'w')
for line in open(args.data):
    info = line.rstrip('\n').split('\t')
    if info[1] in pt_data['participant']:
        sys.stderr.write('Replacing old id ' + info[0] + ' with ' + pt_data['participant'][info[1]] + '\n')
        pt_id_dict[info[0]] = pt_data['participant'][info[1]]
        info[0] = pt_data['participant'][info[1]]
    pt_sheet.write('\t'.join(info) + '\n')
pt_sheet.close()
sys.stderr.write('Finished with patient file.\n')
samp_sheet = open('data_clinical_sample.txt', 'w')

for line in open(args.sample):
    info = line.rstrip('\n').split('\t')
    if info[0] in pt_id_dict:
        sys.stderr.write('Replacing old id ' + info[0] + ' with ' + pt_id_dict[info[0]] + '\n')
        info[0] = pt_id_dict[info[0]]
    samp_sheet.write('\t'.join(info) + '\n')
samp_sheet.close()
sys.stderr.write('Finished with sample file.\n')
