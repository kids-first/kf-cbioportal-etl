import sys
import argparse
import os
import concurrent.futures
import json
import pandas as pd
import numpy as np
from scipy import stats


parser = argparse.ArgumentParser(description='Merge rsem files using cavatica file info.')
parser.add_argument('-d', '--merged_rsem_dir', action='store', dest='merged_rsem_dir', help='merged rsem dir')
parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types '
                                                                                'and data locations')
args = parser.parse_args()

with open(args.config_file) as f:
    config_data = json.load(f)

flist = subprocess.check_output('find ' + args.merged_rsem_dir + ' -name *.rsem_merged.txt', shell=True)
fname_list = flist.decode().split('\n')
if fname_list[-1] == '':
    fname_list.pop()
for fname in fname_list:
    data = pd.read_csv(fname, sep="\t", index_col=0)
    samp_list = list(data.columns)[1:]
    gene_list = list(data.index.values)
    if len(samp_list) > 1:
        new_fname = fname.replace('.rsem_merged.txt', '.rsem_merged_zscore.txt')
        for gene in gene_list:
            data[:gene] = stats.zscore(np.array(data[:gene]), axis=1)
        data.to_csv(new_fname, sep='\t', mode='w', index=False)    
    else:
        sys.stderr.write('Only 1 sample found in ' + fname + ', skipping!\n')
    