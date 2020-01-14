#!/usr/bin/env python3

import sys
import argparse
import concurrent.futures
import re
import json

parser = argparse.ArgumentParser(description='Get all relevant analyzed file outputs from cavatica manifest.')
parser.add_argument('-o', '--output', action='store', dest='output', help='output basename name')
parser.add_argument('-m', '--manifest', action='store', dest='manifest', help='cavatica csv manifest file, may be from step 1a')
parser.add_argument('-b', '--blacklist', action='store', dest='blacklist', help='Optional bs id blacklist')
parser.add_argument('-c', '--config', action='store', dest='config_file', help='json config file with data types and '
                                                                               'data locations')
args = parser.parse_args()
with open(args.config_file) as f:
    config_data = json.load(f)
meta = open(args.manifest)
rel_head = {"Kids First Biospecimen ID Tumor": "tum", "Kids First Biospecimen ID Normal": "norm", "Kids First Biospecimen ID": "rna"}
head = next(meta)
ind = {}
# get index positions of where relevant bs ids would be form tumor dna, normal dna, and rna
header = head.rstrip('\n').split(',')
for i in range(len(header)):
    if header[i] in rel_head:
        ind[rel_head[header[i]]] = i

old_out = open(args.output + '_task_list.txt', 'w')
rna_ext_dict = config_data['rna_ext_list']
rna_ext_list = []
for key in rna_ext_dict:
    rna_ext_list.append(rna_ext_dict[key])
dna_ext_dict = config_data['dna_ext_list']
dna_ext_list = []
for key in dna_ext_dict:
    dna_ext_list.append(dna_ext_dict[key])

skip_dict = {}
if args.blacklist is not None:
    blist = open(args.blacklist)
    for line in blist:
        skip_dict[line.rstrip('\n')] = 0

old_out.write('T/CL BS ID\tNorm BS ID\tTask Name\tTask ID\tAnalyte Type\tRelevant Outputs\tSource Project\n')
cav_dict = {}
for line in meta:
    info = line.rstrip('\n').split(',')
    parts = info[1].split('.')
    ext = ".".join(parts[1:])
    skip = 0
    
    for key in ind:
        if info[ind[key]] in skip_dict:
            skip = 1
            skip_dict[info[ind[key]]] = 1
            sys.stderr.write('BS ID in blacklist.  Skipping ' + line)
    file_id = info[0]
    fname = info[1]
    project = info[2]
    if ext in rna_ext_list and skip == 0:
        bs_id = info[ind['rna']]
        cav_dict[bs_id] = {}
        cav_dict[bs_id]['atype'] = "RNA"
        cav_dict[bs_id]['fname'] = fname
        cav_dict[bs_id]['bs_id'] = bs_id
        cav_dict[bs_id]['project'] = project
    elif ext in dna_ext_list and skip == 0:
        t_bs_id = info[ind['tum']]
        n_bs_id = info[ind['norm']]
        dkey = t_bs_id + n_bs_id
        if dkey not in cav_dict:
            cav_dict[dkey] = {}
            cav_dict[dkey]['atype'] = "DNA"
            cav_dict[dkey]['fname'] = []
            cav_dict[dkey]['t_bs_id'] = t_bs_id
            cav_dict[dkey]['n_bs_id'] = n_bs_id
            cav_dict[dkey]['project'] = project
        cav_dict[dkey]['fname'].append(fname)
for tid in cav_dict:
    if cav_dict[tid]['atype'] == 'RNA':
        old_out.write("\t".join((cav_dict[tid]['bs_id'], "NA","RNA_TASK", tid, "RNA", cav_dict[tid]['fname'], cav_dict[tid]['project'])) + '\n')
    else:
        old_out.write("\t".join((cav_dict[tid]['t_bs_id'], cav_dict[tid]['n_bs_id'],"DNA_TASK", tid, "DNA", ",".join(cav_dict[tid]['fname']), cav_dict[tid]['project'])) + '\n')
old_out.close()