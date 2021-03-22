#!/usr/bin/env python3

import sys
import argparse
import os

parser = argparse.ArgumentParser(description='Quickly add cols to data_clinical sheet')
parser.add_argument('-d', '--datasheet', action='store', dest='data', help='datasheet to add to. Output will be to stdout')
parser.add_argument('-i', '--col_id', action='store', dest='cid', help='Name of column to use for finding ID')
parser.add_argument('-c', '--header', action='store', dest='header', help='tsv with cbio headers to add')
parser.add_argument('-a', '--add', action='store', dest='add', help='csv with BS_ID or PT_ID key and values to add, in same order as header')

args = parser.parse_args()

# tack on new header cols
datasheet = open(args.data)
header = open(args.header)
labels = ""
for head in header:
    entry = next(datasheet)
    sys.stdout.write(entry.rstrip('\n') + "\t" + head)
    labels = entry
header.close()
labels = labels.rstrip('\n').split('\t')
cid = labels.index(args.cid)
# store IDs and vals
value_add = open(args.add)
skip = next(value_add)
add_dict = {}
blanks = ""
for entry in value_add:
    data = entry.rstrip('\n').split(',')
    add_dict[data[0]] = '\t'.join(data[1:])
    # init blanks for when no new value to be added
    if blanks == "":
        for i in range(1, len(data), 1):
            blanks += "\t"
value_add.close()

for line in datasheet:
    info = line.rstrip('\n').split('\t')
    sys.stdout.write(line.rstrip('\n'))
    f = 0
    c = 0
    # IDs might be separeted by ;
    check = info[cid].split(";")
    for i in range(len(check)):
        if check[i] in add_dict:
            f = 1
            c = i
    if f:
        sys.stdout.write("\t" + add_dict[check[c]] + "\n")
    else:
        sys.stdout.write(blanks + "\n")
datasheet.close()
