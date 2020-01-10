
import argparse
import sys
import os
import json


def process_ds(dsname, cav_dict, out):
    try:
        # last item in find will likely be empty after split
        x = 0
        parts = dsname.split('/')
        cbio_proj = parts[-2]
        # project/disease name should be name of directory hosting datasheet
        sys.stderr.write('Processing ' + parts[-2] + ' project' + '\n')
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
                # can't tell if RNA or DNA from datasheet, but DNA will be tum + norm bs id, RNA tum + NA
                id1 = bs_id + "_" + meta[ns_idx]
                id2 = bs_id + "_NA"
                key = id1
                if key in cav_dict:
                    for ftype in cav[id1]:
                        out.write("\t".join([cbio_proj, key, ftype, meta[sa_idx], cav_dict[key][ftype]]) + "\n")
                elif id2 in cav_dict:
                    key = id2
                    for ftype in cav[id1]:
                        out.write("\t".join([cbio_proj, key, ftype, meta[sa_idx], cav_dict[key][ftype]]) + "\n")
                else:
                    sys.stderr.write('WARNING: ' + id1 + ' nor ' id2 + 'found in data sheet entry\n' + entry + 'Check inputs!\n')
    except Exception as e:
        print(e)
        sys.exit()

parser = argparse.ArgumentParser(description='Create temp table with cbio ids, bs ids, and file locations and types for downstream file merging. Should be run one level above datasheets dir.')
parser.add_argument('-c', '--cavatica', action='store', dest='cav',
                    help='file with task info from cavatica (see step 1b)')

args = parser.parse_args()
# config_data dict has most customizable options from json config file
with open(args.config_file) as f:
    config_data = json.load(f)

cav_dict = {}
manifest = open(args.cav)
for line in manifest:
    info = line.rstrip('\n').split(',')
    bs_ids = info[0] + "_" + info[1]
    fnames = info[-2]
    atype = info[4]
    cav_dict[bs_id] = {}
    # will likely have to add fusion ext soon
    for fname in fnames.split(','):
        ftype = 'rsem'
        if atype == 'DNA':
            if fname[-3:] == 'maf':
                ftype = 'maf'
            else:
                ftype =  "cnv"
        cav_dict[bs_id][ftype] = fname

flist = subprocess.check_output('find ./datasheets -name data_clinical_sample.txt', shell=True)

ds_list = flist.decode().split('\n')
if ds_list[-1] == '':
    ds_list.pop()
out = open('cbio_id_fname_table.txt', 'w')
out.write('Cbio project\tT/CL BS ID\tNorm BS ID\tFile Type\tCbio Name\tFile Name\n')
for dpath in ds_list:
    process_ds(dpath, cav_dict, out)
out.close()