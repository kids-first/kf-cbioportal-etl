#!/usr/bin/env python3

import sys
import argparse
import os
import json
import pdb

parser = argparse.ArgumentParser(description='Quickly patch KF PT ID update for PNOC WES')
parser.add_argument('-i', '--input-file', action='store', dest='in_file', help='File to subset')
parser.add_argument('-l', '--id-list', action='store', dest='id_list', help='list of ids seprated by new line')
parser.add_argument('-c', '--columns', action='store', dest='cols', help='csv cols to search')
parser.add_argument('-n', '--num-header-lines', action='store', dest='n', help='number of header lines')
parser.add_argument('-s', '--separator', action='store', dest='s', help='seprator type - \'comma\' or \'tab\'')
parser.add_argument('-o', '--out-dir', action='store', dest='o', help='output dir name - should be different if file being processed is in same dir')
parser.add_argument('-x', '--extra_split', action='store', dest='x', help='split id name on ;, special case')

args = parser.parse_args()
id_file = open(args.id_list)
id_dict = {}
sep = '\t'
if args.s == "comma":
    sep = ','
n = int(args.n)
for line in id_file:
    id_dict[line.rstrip('\n')] = 0

cols = list(map(int, args.cols.split(',')))
#pdb.set_trace()
out_dir = args.o
new_fn = out_dir + os.path.basename(args.in_file)
new_out = open(new_fn, 'w')
in_file = open(args.in_file)
for i in range(n):
    header = next(in_file)
    new_out.write(header)
for line in in_file:
    info = line.rstrip('\n').split(sep)
    if args.x is None:
        for i in cols:
            if info[i] in id_dict:
                id_dict[info[i]] = 1
                new_out.write(line)
    else:
        for i in cols:
            subset = info[i].split(';')
            for j in subset:
                if j in id_dict:
                    id_dict[info[i]] = 1
                    new_out.write(line)
                    break
