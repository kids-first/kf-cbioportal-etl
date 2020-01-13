#!/usr/bin/env python3

import sys
import argparse
import concurrent.futures
import re
import json
import pdb

parser = argparse.ArgumentParser(description='Get all relevant analyzed file outputs from cavatica manifest.')
parser.add_argument('-m', '--manifest', action='store', dest='manifest', help='cavatica manifest file csv list')

args = parser.parse_args()
out_head = 'id,name,project,Kids First Biospecimen ID Tumor,Kids First Biospecimen ID Normal,Kids First Biospecimen ID\n'
out = open('1a_cavatica_merged_manifest.csv', 'w')
out.write(out_head)
for manifests in args.manifest.split(','):
    meta = open(manifests)
    rel_head = ["Kids First Biospecimen ID Tumor", "Kids First Biospecimen ID Normal", "Kids First Biospecimen ID"]
    head = next(meta)
    ind = []

    # get index positions of where relevant bs ids would be form tumor dna, normal dna, and rna
    header = head.rstrip('\n').split(',')
    for i in range(len(rel_head)):
        try:
            ind.append(header.index(rel_head[i]))
        except:
            sys.stderr.write(rel_head[i] + " not in this manifest, assiging 0\n")
            ind.append(0)
    for line in meta:
        info = line.rstrip('\n').split(',')
        out.write(','.join(info[0:3]))
        for idx in ind:
            if idx == 0:
                out.write(',')
            else:
                out.write(',' + info[idx])
        out.write('\n')
out.close()