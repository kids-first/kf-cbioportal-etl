#!/usr/bin/env python3
import sys
import os
import argparse
import json
import subprocess
import concurrent.futures


def process_maf(maf_fn, new_maf, maf_exc, tum_id, norm_id, tid_idx, nid_idx, v_idx, eid_idx, n_idx):
    cur_maf = open(maf_fn)
    next(cur_maf)
    next(cur_maf)
    for maf in cur_maf:
        data = maf.rstrip('\n').split('\t')
        if data[v_idx] not in maf_exc:
            data[tid_idx] = tum_id
            data[nid_idx] = norm_id
            data.pop(eid_idx)
            data.pop((n_idx-1))
            new_maf.write('\t'.join(data) + '\n')
    cur_maf.close()


def process_ds(dsname, c_bs_dict, tid_idx, nid_idx, v_idx, eid_idx, n_idx, print_head):
    try:
        # last item in find will likely be empty after split
        x = 0
        parts = dsname.split('/')
        # project/disease name should be name of directory hosting datasheet
        sys.stderr.write('Processing ' + parts[-2] + ' project' + '\n')
        new_maf = open(out_dir + parts[-2] + '.strelka.vep.filtered.maf', 'w')
        new_maf.write(print_head)
        cur_ds = open(dsname)
        for i in range(0, 4, 1):
            next(cur_ds)
        ds_head = next(cur_ds)
        ds_header = ds_head.rstrip('\n').split('\t')
        sa_idx = ds_header.index('SAMPLE_ID')
        sp_idx = ds_header.index('SPECIMEN_ID')
        ns_idx = ds_header.index('MATCHED_NORMAL_SAMPLE_ID')
        for entry in cur_ds:
            meta = entry.rstrip('\n').split('\t')
            check = meta[sp_idx].split(';')
            for bs_id in check:
                if bs_id in c_bs_dict:
                    sys.stderr.write('Found relevant maf to process for ' + parts[-2] + ' ' + bs_id + '\t' + c_bs_dict[bs_id] + '\t' + meta[sa_idx] + '\n')
                    sys.stderr.flush()
                    process_maf(maf_dir + c_bs_dict[bs_id], new_maf, maf_exc, meta[sa_idx], meta[ns_idx], tid_idx, nid_idx, v_idx, eid_idx, n_idx)
            x += 1
        sys.stderr.write('Completed processing ' + str(x) + ' entries in ' + dsname + '\n')
        new_maf.close()
        cur_ds.close()
    except Exception as e:
        print(e)
        sys.exit()


parser = argparse.ArgumentParser(description='Merge and filter mafs using cavatica task info and datasheets.')
parser.add_argument('-c', '--cavatica', action='store', dest='cav',
                    help='file with task info from cavatica (see step 1)')
parser.add_argument('-i', '--header', action='store', dest='header', help='File with maf header only')
parser.add_argument('-m', '--maf-dir', action='store', dest='maf_dir', help='maf file directory')
parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types and '
                                                                               'data locations')

args = parser.parse_args()
with open(args.config_file) as f:
    config_data = json.load(f)
maf_dir = args.maf_dir
if maf_dir[-1] != '/':
    maf_dir += '/'
cav_fh = open(args.cav)
head = next(cav_fh)
header = head.rstrip('\n').split('\t')
b_idx = header.index('T/CL BS ID')
a_idx = header.index('Analyte Type')
r_idx = header.index('Relevant Outputs')
c_bs_dict = {}
for line in cav_fh:
    info = line.rstrip('\n').split('\t')
    if info[a_idx] == 'DNA':
        files = info[r_idx].split(',')
        for fn in files:
            if fn[-3:] == 'maf':
                c_bs_dict[info[b_idx]] = fn
cav_fh.close()

head_fh = open(args.header)

print_head = next(head_fh)
cur_head = next(head_fh)
cur_header = cur_head.rstrip('\n').split('\t')
eid_idx = cur_header.index('Entrez_Gene_Id')
n_idx = cur_header.index('NCBI_Build')
tid_idx = cur_header.index('Tumor_Sample_Barcode')
nid_idx = cur_header.index('Matched_Norm_Sample_Barcode')
v_idx = cur_header.index('Variant_Classification')
cur_header.pop(eid_idx)
cur_header.pop((n_idx-1))

head_fh.close()
print_head += '\t'.join(cur_header) + '\n'
maf_exc = {"Silent": 0, "Intron": 0, "IGR": 0, "3'UTR": 0, "5'UTR": 0, "3'Flank": 0, "5'Flank": 0}
flist = subprocess.check_output('find ./datasheets -name data_clinical_sample.txt', shell=True)
out_dir = 'merged_mafs/'
try:
    os.mkdir(out_dir)
except:
    sys.stderr.write('output dir already exists\n')

ds_list = flist.decode().split('\n')
if ds_list[-1] == '':
    ds_list.pop()
with concurrent.futures.ProcessPoolExecutor(config_data['cpus']) as executor:
    results = {executor.submit(process_ds, dpath, c_bs_dict, tid_idx, nid_idx, v_idx, eid_idx, n_idx, print_head): dpath for dpath in ds_list}
# for dpath in flist.decode().split('\n'):

sys.stderr.write('Done, check logs\n')
