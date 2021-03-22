#!/usr/bin/env python3
import concurrent.futures
import sys
import pandas as pd
from scipy import stats
import numpy as np

from scripts.file_utils import write_meta_file

def mt_collate_df(config_data,resource):
    sample_id = resource.metadata['sample_id']
    rsem_file = config_data['data_files'] + resource.name    
    current = pd.read_csv(rsem_file, sep="\t", index_col=0)
    cur_subset = current.loc[:, ['FPKM']]
    cur_subset.rename(columns={"FPKM": sample_id}, inplace=True)
    return cur_subset

def add_rsem_file(config_data, sbg_api_client, rsem_resources_by_study):
    sys.stderr.write('Processing for expression files\n')
    seen_list = []
    df_list = []
    for project in rsem_resources_by_study:
        with concurrent.futures.ThreadPoolExecutor(config_data['cpus']) as executor:
            results = {executor.submit(mt_collate_df, config_data, rsem_file): rsem_file for rsem_file in rsem_resources_by_study[project]}
            for result in concurrent.futures.as_completed(results):
                df_list.append(result.result())
    master_tbl = pd.concat(df_list, axis=1)
    del df_list
    # remove duplicate columns/sample_ids
    master_tbl = (master_tbl.loc[:,~master_tbl.columns.duplicated()])
    master_tbl.reset_index(inplace=True)

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

    sample_list = (master_tbl.columns.difference(['Hugo_Symbol']))
    for sym in rpt_sym:
        gene_eval = master_tbl.loc[master_tbl['Hugo_Symbol'] == sym]
        mean_expr = gene_eval[sample_list].mean(axis=1).sort_values(ascending=False)
        to_drop = list(mean_expr.index)[1:]
        master_tbl.drop(to_drop, inplace=True)

    master_tbl.set_index('Hugo_Symbol', inplace=True)
    gene_sym_list = master_tbl.index

    z_scored = stats.zscore(np.log2(np.array(master_tbl + 1)), axis = 1)
    master_zscore_log = pd.DataFrame(z_scored, index=gene_sym_list, columns=sample_list)
    del z_scored
    # this may be memory-intensive for some insane reson...
    sys.stderr.write("Replacing NaN with 0\n")
    sys.stderr.flush()
    master_zscore_log.fillna(0, inplace=True)
    sys.stderr.write('Outputing z scored results\n')
    sys.stderr.flush()

    datasets_dir = config_data['datasets']

    for project in rsem_resources_by_study:
        sampleIds = []
        for resource in rsem_resources_by_study[project]:
            sampleIds.append(resource.metadata['sample_id'])
        expr_fname = '{0}/{1}/data_rna_seq_v2_mrna.txt'.format(datasets_dir, project)
        master_tbl[sampleIds].to_csv(expr_fname, sep='\t', mode='w', index=True)

        zscore_expr_fname = '{0}/{1}/data_rna_seq_v2_mrna_median_Zscores.txt'.format(datasets_dir, project)
        master_zscore_log[sampleIds].to_csv(zscore_expr_fname, sep='\t', mode='w', index=True)

        #Add meta file
        meta_file_path = '{0}/{1}/meta_rna_seq_v2_mrna.txt'.format(datasets_dir,project)
        write_meta_file(project, config_data["metadata"]["rsem"]["counts"]["meta_attr"], meta_file_path)

        meta_file_path = '{0}/{1}/meta_rna_seq_v2_mrna_median_Zscores.txt'.format(datasets_dir, project)
        write_meta_file(project, config_data["metadata"]["rsem"]["zscore"]["meta_attr"], meta_file_path)
