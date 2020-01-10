#!/usr/bin/env python3

import sys
import argparse
import concurrent.futures
import re
import json
import pdb

parser = argparse.ArgumentParser(description='Get all relevant analyzed file outputs from cavatica manifest.')
parser.add_argument('-o', '--output', action='store', dest='output', help='output basename name')
parser.add_argument('-m', '--manifest', action='store', dest='manifest', help='cavatica manifest file csv list')


out_head = 'id,name,project,Kids First Biospecimen ID Tumor,Kids First Biospecimen ID Normal,Kids First Biospecimen ID\n'
out = open('1a_cavatica_merged_manifest.csv', 'w')
for manifests in arg.manifest.split(','):
    meta = open(manifests)
    rel_head = {"Kids First Biospecimen ID Tumor": "tum", "Kids First Biospecimen ID Normal": "norm", "Kids First Biospecimen ID": "rna"}
    head = next(meta)
    ind = []
    # get index positions of where relevant bs ids would be form tumor dna, normal dna, and rna
    header = head.rstrip('\n').split(',')
    for i in range(len(header)):
        if header[i] in rel_head:
            ind.append[i]
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