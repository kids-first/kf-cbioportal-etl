#!/usr/bin/env python3

import sys
import argparse
import os
import concurrent.futures
import gzip
import pandas as pd
import numpy as np
from scipy import stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge rsem files using cavatica file info.')
    parser.add_argument('-f', '--file', action='store', dest='in_file',
                    help='log2 transformed FPKM values')
    parser.add_argument('-r', '--rsem-dir', action='store', dest='rsem_dir', help='rsem file directory')
    args = parser.parse_args()


    rsem_dir = args.rsem_dir

    if rsem_dir[-1] != '/':
        rsem_dir += '/'
    
    master_tbl = pd.read_csv(args.in_file, sep="\t")
    
    master_tbl.set_index('Hugo_Symbol', inplace=True)
    gene_sym_list = master_tbl.index
    sample_list=list(master_tbl.columns.values)
    sys.stderr.write('Calculating z scores\n')
    sys.stderr.flush()

    z_scored = stats.zscore(np.log2(np.array(master_tbl + 1)), axis = 1)
    del master_tbl
    master_zscore_log = pd.DataFrame(z_scored, index=gene_sym_list, columns=sample_list)
    master_zscore_log.fillna(0, inplace=True)
    master_zscore_log.to_csv("cptac_published.zscore.log2_fpkm.tsv", sep='\t', mode='w', index=True)
