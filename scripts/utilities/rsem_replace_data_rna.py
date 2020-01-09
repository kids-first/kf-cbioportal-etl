#!/usr/bin/env python3

import sys
import argparse
import os

parser = argparse.ArgumentParser(description='Updates outdated data in data_rna_seq sheet from new rsem_merged file')
parser.add_argument('-m', '--mapFileLoc', action='store', dest='map', help='RNA sample name map file')
parser.add_argument('-r', '--rsem', action='store', dest='rsem', help='rsem file with new data')
parser.add_argument('-d', '--data', action='store', dest='data', help='data_rna_seq_v2_mrna.txt sheet to update')
parser.add_argument('-o', '--outdir', action='store', dest='outdir', help='outputdir')

args = parser.parse_args()
map_dict = {}
for line in open(args.map):
    (new_name, old_name) = line.rstrip('\n').split('\t')
    map_dict[old_name] = new_name

rsem_data = {}

rsem_file = open(args.rsem)
head = next(rsem_file)
slist = []
header = head.rstrip('\n').split('\t')
for line in rsem_file:
    info = line.rstrip('\n').split('\t')
    rsem_data[info[1]] = {}
    for i in range(2, len(info)):
        samp = map_dict[header[i]]
        slist.append(samp)
        rsem_data[info[1]][samp] = info[i]

rsem_file.close()

outdir = args.outdir + '/' + os.path.dirname(args.data)
if not os.path.isdir(outdir):
    os.mkdir(outdir)
update = open(outdir + '/data_rna_seq_v2_mrna.txt', 'w')

old_data = open(args.data)
head = next(old_data)
update.write(head)
header = head.rstrip('\n').split('\t')

h_i = []
for i in range(1, len(header)):
    if header[i] in slist:
        h_i.append(i)
        sys.stderr.write('Found sample ' + header[i] + ' to update at index ' + str(i) + '\n')
for line in old_data:
    info = line.rstrip('\n').split('\t')
    for i in h_i:
        info[i] = rsem_data[info[0]][header[i]]
    update.write('\t'.join(info) + '\n')
update.close()
