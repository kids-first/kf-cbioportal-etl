#!/usr/bin/env python3

import sys
import argparse
import os
import concurrent.futures
import gzip
import json
import pandas as pd
import numpy as np
from scipy import stats


def mt_process_master(rsem_file):
    try:
        current = pd.read_csv(rsem_dir + rsem_list[i], sep="\t", index_col=0)
        sample = rna_subset.loc[rna_subset['File_Name'] == rsem_list[i], 'Cbio_Tumor_Name'].iloc[0]
        sample_list.append(sample)
        master_tbl = pd.concat([master_tbl, current[['FPKM']]], axis=1)
    except Exception as e:
        sys.stderr.write('Failed to concat ' + rsem_file + '\n')
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge rsem files using cavatica file info.')
    parser.add_argument('-t', '--table', action='store', dest='table',
                    help='Table with cbio project, kf bs ids, cbio IDs, and file names')
    parser.add_argument('-r', '--rsem-dir', action='store', dest='rsem_dir', help='rsem file directory')
    parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types '
                                                                                   'and data locations')
    args = parser.parse_args()

    with open(args.config_file) as f:
        config_data = json.load(f)
    rsem_dir = args.rsem_dir

    if rsem_dir[-1] != '/':
        rsem_dir += '/'

    # file_meta_dict = get_file_metadata(args.table, 'rsem')
    out_dir = 'merged_rsem/'
    try:
        os.mkdir(out_dir)
    except:
        sys.stderr.write('output dir already exists\n')
        sys.stderr.flush()
    
    # Use cbio table to create master table to choose represntative gene symbol, calc z scores, output both
    all_file_meta = pd.read_csv(args.table, sep="\t")
    rna_subset = file_meta.loc[file_meta['File_Type'] == 'rsem']
    rsem_list = rns_subset['File_Name'].to_list()
    init_tbl = pd.read_csv(rsem_dir + rsem_list[0], sep="\t", index_col=0)
    # Get cbio name to rename columns in mster table
    sample_list = []
    sample = rna_subset.loc[rna_subset['File_Name'] == rsem_list[0], 'Cbio_Tumor_Name'].iloc[0]
    sample_list.append(sample)
    master_tbl = init_tbl[['FPKM']]
    # master_tbl.rename(columns={"FPKM": sample})
    sys.stderr.write('Creating merged rsem table\n')
    x = 1
    m = 100

    with concurrent.futures.ThreadPoolExecutor(config_data['threads']) as executor:
        results = {executor.submit(mt_process_master, rsem_list[i],): rsem_list[i] for i in range(1, len(rsem_list), 1)}
        for result in concurrent.futures.as_completed(results):
            if x % m == 0:
                sys.stderr.write('Merged ' + str(x) + ' files\n')
                sys.stderr.flush()
    master_tbl.columns = sample_list
    
    sys.stderr.write('Merge complete, picking representative gene transcripts\n')
    master_tbl.reset_index(inplace=True)
    # drop ens_id_sym and replace with just sym
    gene_id_sym_list = master_tbl['gene_id'].to_list()
    gene_sym_list = []
    for entry in gene_id_sym_list:
        parts = entry.split('_')
        gene_sym = '_'.join(parts[1:])
        gene_sym_list.append(gene_sym)
    master_tbl['Hugo_Symbol'] = gene_sym_list
    master_tbl.drop(columns =["gene_id"], inplace = True)
    dup_hugo_tbl = master_tbl[master_tbl.duplicated(['Hugo_Symbol'])]
    rpt_sym = dup_hugo_tbl.Hugo_Symbol.unique()
    for sym in rpt_sym:
        gene_eval = master_tbl.loc[master_tbl['Hugo_Symbol'] == sym]
        mean_expr = gene_eval[sample_list].mean(axis=1).sort_values(ascending=False)
        to_dump = list(mean_expr.index)[1:]
        master_tbl.drop(to_dump, inplace=True)
    
    sys.stderr.write('Calculating z scores\n')
    master_tbl.set_index('Hugo_Symbol', inplace=True)
    gene_sym_list = master_tbl.index
    log_tf = np.log(np.array(master_tbl))
    z_scored = stats.zscore(log_tf, axis = 1)
    master_zscore_log = pd.DataFrame(z_scored, index=gene_sym_list, columns=sample_list)
    master_zscore_log.fillna('NA', inplace=True)
    
    sys.stderr.write('Outputing project files\n')
    project_list = rna_subset.Cbio_project.unique()
    for project in project_list:
        sub_sample_list = list(rna_subset.loc[rna_subset['Cbio_project'] == project_list[0], 'Cbio_Tumor_Name'])
        expr_fname = out_dir + project + '.rsem_merged.txt'
        master_tbl[sub_sample_list].to_csv(expr_fname, sep='\t', mode='w', index=True)
        zscore_fname = out_dir + project + '.rsem_merged_zscore.txt'
        master_zscore_log[sub_sample_list].to_csv(zscore_fname, sep='\t', mode='w', index=True)