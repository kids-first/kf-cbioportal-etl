#!/usr/bin/env python3

import sevenbridges as sbg
import sys
import argparse
import concurrent.futures
from sevenbridges.http.error_handlers import rate_limit_sleeper, maintenance_sleeper
import re
import json
import pdb

parser = argparse.ArgumentParser(description='Get all relevant analyzed file outputs from cavatica manifest.')
parser.add_argument('-o', '--output', action='store', dest='output', help='output basename name')
parser.add_argument('-m', '--manifest', action='store', dest='manifest', help='cavatica csv manifest file')
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

mafs = open(args.output + '_maf_file_list.txt', 'w')
cnvs = open(args.output + '_cnv_file_list.txt', 'w')
rsem = open(args.output + '_rsem_file_list.txt', 'w')
old_out = open(args.output + '_task_list.txt', 'w')
rna_ext_list = config_data['rna_ext_list']
dna_ext_list = config_data['dna_ext_list']

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

    if ext in rna_ext_list and skip == 0:
        bs_id = info[ind['rna']]
        rsem.write("\t".join((info[0], info[1], bs_id)) + "\n")
        cav_dict[bs_id] = {}
        cav_dict[bs_id]['atype'] = "RNA"
        cav_dict[bs_id]['fname'] = info[1]
        cav_dict[bs_id]['bs_id'] = bs_id
    elif ext in dna_ext_list and skip == 0:
        t_bs_id = info[ind['tum']]
        n_bs_id = info[ind['norm']]
        dkey = t_bs_id + n_bs_id
        if re.search('.maf', info[1]):
            mafs.write("\t".join((info[0], info[1])) + "\n")
        else:
            cnvs.write("\t".join((info[0], info[1])) + "\n")
        if dkey not in cav_dict:
            cav_dict[dkey] = {}
            cav_dict[dkey]['atype'] = "DNA"
            cav_dict[dkey]['fname'] = []
            cav_dict[dkey]['t_bs_id'] = t_bs_id
            cav_dict[dkey]['n_bs_id'] = n_bs_id
        cav_dict[dkey]['fname'].append(info[1])
mafs.close()
cnvs.close()
rsem.close()
proj = config_data['cancStudyID']
for tid in cav_dict:
    if cav_dict[tid]['atype'] == 'RNA':
        old_out.write("\t".join((cav_dict[tid]['bs_id'], "NA","RNA_TASK", tid, "RNA", cav_dict[tid]['fname'], proj)) + '\n')
    else:
        old_out.write("\t".join((cav_dict[tid]['t_bs_id'], cav_dict[tid]['n_bs_id'],"DNA_TASK", tid, "DNA", ",".join(cav_dict[tid]['fname']), proj)) + '\n')
old_out.close()