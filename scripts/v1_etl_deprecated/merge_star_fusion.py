#!/usr/bin/env python3
import sys
import os
import argparse
import json

parser = argparse.ArgumentParser(description='Merge and add sample id to STAR Fusion output.')
parser.add_argument('-d', '--dir', action='store', dest='dir', help='star fusion file directory')
parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types and '
                                                                               'data locations')
parser.add_argument('-i', '--header', action='store', dest='header', help='File with fusion header only')
parser.add_argument('-m', '--manifest', action='store', dest='manifest', help='star fusion cavatica manifest file')

args = parser.parse_args()
with open(args.config_file) as f:
    config_data = json.load(f)
h_file = open(args.header)
header = next(h_file)
header = "tumor_id\t" + header
merged = open('merged_fusions.tsv', 'w')
merged.write(header)
map_dict = {}
for line in open(config_data['mapFileLoc']):
    (samp_id, bs_id) = line.rstrip('\n').split('\t')
    map_dict[bs_id] = samp_id



manifest = open(args.manifest)
head = next(manifest)
header = head.rstrip('\n').split(',')
b_idx = header.index('Kids First Biospecimen ID')

for line in manifest:
    info = line.rstrip('\n').split(',')
    bs_id = info[b_idx]
    if bs_id in map_dict:
        samp_id = map_dict[bs_id]
        floc = args.dir + "/" + info[1]
        fus = open(floc)
        skip = next(fus)
        for fdata in fus:
            merged.write(samp_id + "\t" + fdata)
merged.close()