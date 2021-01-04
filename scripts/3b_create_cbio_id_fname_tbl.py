#!/usr/bin/env python3
import argparse
import sys
import os
import json
import subprocess


def process_ds(dsname, cav_dict, out):
    try:
        # last item in find will likely be empty after split
        x = 0
        parts = dsname.split('/')
        cbio_proj = parts[-2]
        # project/disease name should be name of directory hosting datasheet
        sys.stderr.write('Processing ' + cbio_proj + ' project' + '\n')
        cur_ds = open(dsname)
        for i in range(0, 4, 1):
            next(cur_ds)
        ds_head = next(cur_ds)
        ds_header = ds_head.rstrip('\n').split('\t')
        sa_idx = ds_header.index('SAMPLE_ID')
        sp_idx = ds_header.index('SPECIMEN_ID')
        ns_idx = ds_header.index('MATCHED_NORMAL_SPECIMEN_ID')
        for entry in cur_ds:
            meta = entry.rstrip('\n').split('\t')
            check = meta[sp_idx].split(';')
            for bs_id in check:
                # can't tell if RNA or DNA from datasheet, but DNA will be tum + norm bs id, RNA tum + NA
                id1 = bs_id + "\t" + meta[ns_idx]
                id2 = bs_id + "\tNA"
                key = id1
                if key in cav_dict:
                    for ftype in cav_dict[key]:
                        out.write("\t".join([cbio_proj, key, ftype, meta[sa_idx], meta[ns_idx], cav_dict[key][ftype]]) + "\n")
                elif id2 in cav_dict:
                    key = id2
                    for ftype in cav_dict[key]:
                        out.write("\t".join([cbio_proj, key, ftype, meta[sa_idx], 'NA', cav_dict[key][ftype]]) + "\n")
                else:
                    sys.stderr.write('WARNING: ' + id1 + ' nor ' + id2 + ' found in data sheet entry\n' + entry + 'Check inputs!\n')
    except Exception as e:
        print(e)
        sys.exit()

parser = argparse.ArgumentParser(description='Create temp table with cbio ids, bs ids, and file locations and types for downstream file merging. Should be run one level above datasheets dir.')
parser.add_argument('-c', '--cavatica', action='store', dest='cav',
                    help='file with task info from cavatica (see step 1b)')

args = parser.parse_args()

cav_dict = {}
manifest = open(args.cav)
head = next(manifest)
for line in manifest:
    info = line.rstrip('\n').split('\t')
    bs_ids = info[0] + "\t" + info[1]
    fnames = info[-2]
    atype = info[4]
    cav_dict[bs_ids] = {}
    # will likely have to add fusion ext soon
    for fname in fnames.split(','):
        ftype = 'rsem'
        if atype == 'DNA':
            if fname[-3:] == 'maf':
                ftype = 'maf'
            else:
                ftype =  "cnv"
        elif fname[-3:] == 'tsv':
            ftype = 'fusion'
        cav_dict[bs_ids][ftype] = fname

flist = subprocess.check_output('find /Users/kalletlak/Documents/datasheets -name data_clinical_sample.txt', shell=True)

ds_list = flist.decode().split('\n')
if ds_list[-1] == '':
    ds_list.pop()
out = open('cbio_id_fname_table.txt', 'w')
out.write('Cbio_project\tT_CL_BS_ID\tNorm_BS_ID\tFile_Type\tCbio_Tumor_Name\tCbio_Matched_Normal_Name\tFile_Name\n')
for dpath in ds_list:
    process_ds(dpath, cav_dict, out)
out.close()