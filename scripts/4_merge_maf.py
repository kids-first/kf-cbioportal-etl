#!/usr/bin/env python3
import sys
import os
import argparse
import json
import subprocess
from get_file_metadata_helper import get_file_metadata
import concurrent.futures


def process_maf(maf_fn, new_maf, maf_exc, tum_id, norm_id, tid_idx, h_idx, nid_idx, v_idx, eid_idx):
    cur_maf = open(maf_fn)
    next(cur_maf)
    next(cur_maf)
    for maf in cur_maf:
        data = maf.rstrip('\n').split('\t')
        # Want to allow TERT promoter as exception to exlcusion rules
        if data[v_idx] not in maf_exc or (data[h_idx] == "TERT" and data[v_idx] == "5'Flank"):
            data[tid_idx] = tum_id
            data[nid_idx] = norm_id
            data.pop(eid_idx)
            new_maf.write('\t'.join(data) + '\n')
    cur_maf.close()


def process_tbl(cbio_dx, file_meta_dict, tid_idx, h_idx, nid_idx, v_idx, eid_idx, print_head):
    try:
        x = 0
        # project/disease name should be name of directory hosting datasheet
        sys.stderr.write('Processing ' + cbio_dx + ' project' + '\n')
        new_maf = open(out_dir + cbio_dx + ".maf", 'w')
        new_maf.write(print_head)
        for cbio_tum_id in file_meta_dict[cbio_dx]:
            cbio_norm_id = file_meta_dict[cbio_dx][cbio_tum_id]['cbio_norm_id']
            fname = file_meta_dict[cbio_dx][cbio_tum_id]['fname']
            sys.stderr.write('Found relevant maf to process for ' + ' ' + cbio_tum_id + ' ' + cbio_norm_id + ' '
            + file_meta_dict[cbio_dx][cbio_tum_id]['kf_tum_id'] + ' ' + file_meta_dict[cbio_dx][cbio_tum_id]['kf_norm_id'] + ' ' + fname + '\n')
            sys.stderr.flush()
            process_maf(maf_dir + fname, new_maf, maf_exc, cbio_tum_id, cbio_norm_id, tid_idx, h_idx, nid_idx, v_idx, eid_idx)
            x += 1
        sys.stderr.write('Completed processing ' + str(x) + ' entries in ' + cbio_dx + '\n')
        new_maf.close()
    except Exception as e:
        print(e)
        sys.exit()


parser = argparse.ArgumentParser(description='Merge and filter mafs using cavatica task info and datasheets.')
parser.add_argument('-t', '--table', action='store', dest='table',
                    help='Table with cbio project, kf bs ids, cbio IDs, and file names')
parser.add_argument('-i', '--header', action='store', dest='header', help='File with maf header only')
parser.add_argument('-m', '--maf-dir', action='store', dest='maf_dir', help='maf file directory')
parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types and '
                                                                               'data locations')

args = parser.parse_args()
with open(args.config_file) as f:
    config_data = json.load(f)
# get maf file ext
maf_dir = args.maf_dir
if maf_dir[-1] != '/':
    maf_dir += '/'
file_meta_dict = get_file_metadata(args.table, 'maf')
head_fh = open(args.header)

print_head = next(head_fh)
cur_head = next(head_fh)
cur_header = cur_head.rstrip('\n').split('\t')
eid_idx = cur_header.index('Entrez_Gene_Id')
tid_idx = cur_header.index('Tumor_Sample_Barcode')
nid_idx = cur_header.index('Matched_Norm_Sample_Barcode')
v_idx = cur_header.index('Variant_Classification')
h_idx = cur_header.index('Hugo_Symbol')
cur_header.pop(eid_idx)

head_fh.close()
print_head += '\t'.join(cur_header) + '\n'
maf_exc = {"Silent": 0, "Intron": 0, "IGR": 0, "3'UTR": 0, "5'UTR": 0, "3'Flank": 0, "5'Flank": 0, "RNA": 0}
out_dir = 'merged_mafs/'
try:
    os.mkdir(out_dir)
except:
    sys.stderr.write('output dir already exists\n')

with concurrent.futures.ProcessPoolExecutor(config_data['cpus']) as executor:
    results = {executor.submit(process_tbl, cbio_dx, file_meta_dict, tid_idx, h_idx, nid_idx, v_idx, eid_idx, print_head): cbio_dx for cbio_dx in file_meta_dict}

sys.stderr.write('Done, check logs\n')