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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge rsem files using cavatica file info.')
    parser.add_argument('-t', '--table', action='store', dest='table',
                    help='Table with cbio project, kf bs ids, cbio IDs, and file names')
    parser.add_argument('-r', '--rsem-dir', action='store', dest='rsem_dir', help='rsem file directory')
    parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types '
                                                                                   'and data locations')
    args = parser.parse_args()

    def mt_collate_df(rsem_file):
        try:
            sample = rna_subset.loc[rna_subset['File_Name'] == rsem_file, 'Cbio_Tumor_Name'].iloc[0]
            if sample not in seen_list:
                sample_list.append(sample)
                current = pd.read_csv(rsem_dir + rsem_file, sep="\t", index_col=0)
                cur_subset = current[['FPKM']]
                cur_subset.rename(columns={"FPKM": sample}, inplace=True)
                df_list.append(cur_subset)
                seen_list.append(sample)
            else:
                sys.stderr.write(sample + " was already seen, likely due to mulitple dx, skipping!\n")
            return 1
        except Exception as e:
            sys.stderr.write(str(e) + "\nFailed to process " + rsem_file)

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
    rna_subset = all_file_meta.loc[all_file_meta['File_Type'] == 'rsem']
    rsem_list = rna_subset['File_Name'].to_list()
    # init_tbl = pd.read_csv(rsem_dir + rsem_list[0], sep="\t", index_col=0)
    # Get cbio name to rename columns in master table
    sample_list = []
    # master_tbl = init_tbl[['FPKM']]
    # master_tbl.rename(columns={"FPKM": sample}, inplace=True)
    sys.stderr.write('Creating merged rsem table\n')
    x = 1
    m = 50

    df_list = []
    # some samples have more than one dx, need only process once
    seen_list = []
    with concurrent.futures.ThreadPoolExecutor(16) as executor:
        results = {executor.submit(mt_collate_df, rsem_list[i]): rsem_list[i] for i in range(len(rsem_list))}
        for result in concurrent.futures.as_completed(results):
            if x % m == 0:
                sys.stderr.write('Merged ' + str(x) + ' files\n')
                sys.stderr.flush()
            x += 1
    master_tbl = pd.concat(df_list, axis=1)
    del df_list

    sys.stderr.write('Merge complete, picking representative gene transcripts\n')
    sys.stderr.flush()
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
    
    master_tbl.set_index('Hugo_Symbol', inplace=True)
    gene_sym_list = master_tbl.index
    project_list = rna_subset.Cbio_project.unique()
    
    sys.stderr.write('Outputting FPKM expression results\n')
    sys.stderr.flush()
    for project in project_list:
        sub_sample_list = list(rna_subset.loc[rna_subset['Cbio_project'] == project, 'Cbio_Tumor_Name'])
        expr_fname = out_dir + project + '.rsem_merged.txt'
        master_tbl[sub_sample_list].to_csv(expr_fname, sep='\t', mode='w', index=True)
    
    sys.stderr.write('Calculating z scores\n')
    sys.stderr.flush()

    z_scored = stats.zscore(np.log2(np.array(master_tbl + 1)), axis = 1)
    del master_tbl
    master_zscore_log = pd.DataFrame(z_scored, index=gene_sym_list, columns=sample_list)
    # this may be memory-intensive for some insane reson...
    sys.stderr.write("Replacing NaN with 0\n")
    sys.stderr.flush()
    master_zscore_log.fillna(0, inplace=True)
    sys.stderr.write('Outputing z scored results\n')
    sys.stderr.flush()
    
    for project in project_list:
        sub_sample_list = list(rna_subset.loc[rna_subset['Cbio_project'] == project, 'Cbio_Tumor_Name'])
        zscore_fname = out_dir + project + '.rsem_merged_zscore.txt'
        master_zscore_log[sub_sample_list].to_csv(zscore_fname, sep='\t', mode='w', index=True)