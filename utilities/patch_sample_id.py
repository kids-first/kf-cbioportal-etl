#!/usr/bin/env python3

import sys
import argparse
import os

parser = argparse.ArgumentParser(description='Quickly patch sample id from step2 with data warehouse entries')
parser.add_argument('-i', '--info', action='store', dest='info', help='A sample info sheet from step 2')
parser.add_argument('-d', '--data', action='store', dest='data', help='Data Warehouse tsv with biospecimen ID, master aliquot id, and sample id')

new_samp = {}
args = parser.parse_args()

data_in = open(args.data)
head = next(data_in)
header = head.rstrip('\n').split('\t')
b_idx = header.index('kf_biospecimen_id')
s_idx = header.index('sample_id')
a_idx = header.index('master_parent_aliquot_id')

for line in data_in:
    meta = line.rstrip('\n').split('\t')
    new_samp[meta[b_idx]] = meta[s_idx] + "_" + meta[a_idx]
data_in.close()

info_in = open(args.info)
head = next(info_in)
header = head.rstrip('\n').split('\t')
b_idx = header.index('BS_ID')
s_idx = header.index('external_sample_id')
a_idx = header.index('analyte_type')

sys.stdout.write(head)
for line in info_in:
    meta = line.rstrip('\n').split('\t')
    if meta[b_idx] in new_samp:
        meta[s_idx] = new_samp[meta[b_idx]]
    else:
        sys.stderr.write("WARN: Entry for " + meta[b_idx] + " " + meta[s_idx] + " not found, using original\n")
    if meta[a_idx] == "Not Applicable":
        meta[a_idx] = "DNA"
    print("\t".join(meta))
info_in.close()

