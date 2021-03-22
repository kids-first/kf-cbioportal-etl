#!/usr/bin/env python3

import sys
import argparse

parser = argparse.ArgumentParser(description='Quickly add cols to data_clinical sheet')
parser.add_argument('-c', '--col-names-file', action='store', dest='col_file', help='File the new line seprated header names to null out')
parser.add_argument('-m', '--maf_in', action='store', dest='maf_in', help='MAF file with cols to null. Output will be to stdout')

args = parser.parse_args()

col_file = args.col_file
maf = args.maf_in

col_names = []
for line in open(col_file):
    col_names.append(line.rstrip('\n'))

col_idx = []
maf_rd = open(maf)
head = next(maf_rd)
sys.stdout.write(head)
head = next(maf_rd)
header = head.rstrip('\n').split('\t')
for col in col_names:
    col_idx.append(header.index(col))
sys.stdout.write(head)
for data in maf_rd:
    entry = data.rstrip('\n').split('\t')
    for i in col_idx:
        entry[i] = ""
    sys.stdout.write("\t".join(entry) + "\n")
