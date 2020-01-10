#!/usr/bin/env python3

import sys
import argparse
import os
import concurrent.futures
import gzip
import json


def process_rsem_file(rsem_info, rsem_dir, fpkm_cts):
    try:
        rsem_path = rsem_dir + rsem_info[1]
        bs_id = rsem_info[2]
        if bs_id in mapping_dict:
            samp_id = mapping_dict[bs_id]
            samp_list.append(samp_id)
            cur = gzip.open(rsem_dir + rsem_info[1])
            next(cur)
            for line in cur:
                data = line.decode().rstrip('\n').split('\t')
                gene_info = data[0].split('_')
                ensg_id = gene_info[0]
                fpkm_cts[ensg_id][samp_id] = data[-1]
            cur.close()
        else:
            sys.stderr.write(rsem_path + ' not in mapping file, skipping!\n')
            sys.stderr.flush()
    except Exception as e:
        sys.stderr.write('Error while processing ' + rsem_path + ', check that the file exists\n')
        print(e)
        sys.exit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge rsem files using cavatica file info.')
    parser.add_argument('-f', '--rsem_file_list', action='store', dest='flist',
                        help='File from step 1b, rsem file ids and file names')
    parser.add_argument('-r', '--rsem-dir', action='store', dest='rsem_dir', help='rsem file directory')
    parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types '
                                                                                   'and data locations')
    args = parser.parse_args()

    with open(args.config_file) as f:
        config_data = json.load(f)
    id_list = []
    id_file = open(config_data['ens_gene_list'])
    next(id_file)
    gene_dict = {}
    for line in id_file:
        (ens, sym) = line.rstrip('\n').split('\t')
        id_list.append(ens)
        gene_dict[ens] = sym

    mapping_dict = {}
    samp_list = []
    for line in open(config_data['mapFileLoc']):
        info = line.rstrip('\n').split('\t')
        mapping_dict[info[1]] = info[0]
    rsem_fh = open(args.flist)
    rsem_dir = args.rsem_dir

    if rsem_dir[-1] != '/':
        rsem_dir += '/'
    rsem_list =[]
    for line in rsem_fh:
        rsem_list.append(line.rstrip('\n').split('\t'))
    rsem_fh.close()
    

    out_dir = 'merged_rsem/'
    try:
        os.mkdir(out_dir)
    except:
        sys.stderr.write('output dir already exists\n')
        sys.stderr.flush()
    fpkm_cts = {}
    for ens_id in id_list:
        fpkm_cts[ens_id] = {}
    # working with shared dict with "global" scope, no need to pass as param or may be overwritten

    with concurrent.futures.ThreadPoolExecutor(config_data['threads']) as executor:
        i = 1
        n = 50
        merge_results = {executor.submit(process_rsem_file, rsem_info, rsem_dir, fpkm_cts): rsem_info for rsem_info in rsem_list}
        for processed in concurrent.futures.as_completed(merge_results):
            if i % n == 0:
                sys.stderr.write('Processed ' + str(i) + ' input files\n')
            i += 1
    samp_list.sort()
    rsem_merge = open(out_dir + 'rsem_merged.txt', 'w')
    rsem_merge.write('gene_id' + '\t' + 'gene_symbol' + '\t' + '\t'.join(samp_list) + '\n')
    for ens_id in id_list:
        rsem_merge.write(ens_id + '\t' + gene_dict[ens_id])
        for samp_id in samp_list:
            rsem_merge.write('\t' + fpkm_cts[ens_id][samp_id])
        rsem_merge.write('\n')
    rsem_merge.close()

    sys.stderr.write('Done, check logs\n')
