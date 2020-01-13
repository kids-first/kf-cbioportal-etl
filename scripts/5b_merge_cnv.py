#!/usr/bin/env python3
import sys
import os
import argparse
import json
import subprocess
import concurrent.futures


def process_cnv(cnv_fn, cur_cnv_dict, samp_id):
    for entry in open(cnv_fn):
        (gene, gid, value) = entry.rstrip('\n').split('\t')
        if gene not in cur_cnv_dict:
            cur_cnv_dict[gene] = {}
        if samp_id in cur_cnv_dict[gene]:
            sys.stderr.write('ERROR: Sample ID ' + samp_id + ' already processed.  Back to the drawing board!')
            exit(1)
        else:
            cur_cnv_dict[gene][samp_id] = value
    return cur_cnv_dict


def process_ds(dpath, c_bs_dict):
    try:
        if dpath != '':
            parts = dpath.split('/')
            # project/disease name should be name of directory hosting datasheet
            sys.stderr.write('Processing ' + parts[-2] + ' project' + '\n')
            new_cnv = open(out_dir + parts[-2] + '.predicted_cnv.txt', 'w')
            cur_cnv_dict = {}
            s_list = []
            cur_ds = open(dpath)
            for i in range(0, 4, 1):
                next(cur_ds)
            ds_head = next(cur_ds)
            ds_header = ds_head.rstrip('\n').split('\t')
            sa_idx = ds_header.index('SAMPLE_ID')
            sp_idx = ds_header.index('SPECIMEN_ID')
            for line in cur_ds:
                info = line.rstrip('\n').split('\t')
                check = info[sp_idx].split(';')
                for bs_id in check:
                    if bs_id in c_bs_dict:
                        sys.stderr.write('Found relevant cnv to process ' + bs_id + '\t' + c_bs_dict[bs_id] + '\t' + info[sa_idx] + '\n')
                        sys.stderr.flush()
                        s_list.append(info[sa_idx])
                        cur_cnv_dict = process_cnv(cnv_dir + c_bs_dict[bs_id] + suffix, cur_cnv_dict, info[sa_idx])
            new_cnv.write('Hugo_Symbol\t' + '\t'.join(s_list) + '\n')
            for gene in cur_cnv_dict:
                new_cnv.write(gene)
                for samp in s_list:
                    if samp in cur_cnv_dict[gene]:
                        new_cnv.write('\t' + cur_cnv_dict[gene][samp])
                    else:
                        new_cnv.write('\t2')
                new_cnv.write('\n')
            new_cnv.close()
            cur_ds.close()
    except Exception as e:
        print(e)
        sys.exit(1)


parser = argparse.ArgumentParser(description='Merge cnv files using cavatica task info and datasheets.')
parser.add_argument('-c', '--cavatica', action='store', dest='cav',
                    help='file with task info from cavatica (see step 1)')
parser.add_argument('-n', '--cnv-dir-gene', action='store', dest='cnv_dir', help='cnv as gene file directory')
parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types and '
                                                                               'data locations')

args = parser.parse_args()
with open(args.config_file) as f:
    config_data = json.load(f)
cav_fh = open(args.cav)
cnv_dir = args.cnv_dir
if cnv_dir[-1] != '/':
    cnv_dir += '/'

suffix = '.CNVs.Genes.copy_number'

flist = subprocess.check_output('find ./datasheets -name data_clinical_sample.txt', shell=True)
out_dir = 'merged_cnvs/'
try:
    os.mkdir(out_dir)
except:
    sys.stderr.write('output dir already exists\n')

with concurrent.futures.ProcessPoolExecutor(config_data['cpus']) as executor:
    results = {executor.submit(process_ds, dpath, c_bs_dict): dpath for dpath in flist.decode().split('\n')}

sys.stderr.write('Done, check logs\n')
