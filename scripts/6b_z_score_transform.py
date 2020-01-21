import sys
import argparse
import os
import concurrent.futures
import json
import subprocess
import pandas as pd
import numpy as np
from scipy import stats


parser = argparse.ArgumentParser(description='Merge rsem files using cavatica file info.')
parser.add_argument('-d', '--merged_rsem_dir', action='store', dest='merged_rsem_dir', help='merged rsem dir')
parser.add_argument('-t', '--table', action='store', dest='table',
                    help='Table with cbio project, kf bs ids, cbio IDs, and file names, ' +
                    'give only if you want to normalize first across all tables, then split')
parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types '
                                                                                'and data locations')
args = parser.parse_args()

with open(args.config_file) as f:
    config_data = json.load(f)

flist = subprocess.check_output('find ' + args.merged_rsem_dir + ' -name *.rsem_merged.txt', shell=True)
fname_list = flist.decode().split('\n')
if fname_list[-1] == '':
    fname_list.pop()

if args.table == None:
    for fname in fname_list:
        data = pd.read_csv(fname, sep="\t", index_col=0)
        header = list(data.columns)
        gene_list = list(data.index.values)
        if len(header) > 2:
            new_fname = fname.replace('.rsem_merged.txt', '.rsem_merged_zscore.txt')
            z_scored = stats.zscore(np.array(data), axis=1)
            new_data = pd.DataFrame(z_scored, index=gene_list, columns=header)
            new_data.fillna('NA', inplace=True)
            sys.stderr.write('Outputting results to ' + new_fname + '\n')
            new_data.to_csv(new_fname, sep='\t', mode='w', index=True)
        else:
            sys.stderr.write('Only 1 sample found in ' + fname + ', skipping!\n')
else:
    mega_tbl = pd.read_csv(fname_list[0], sep="\t", index_col=0)
    mega_header = list(mega_tbl.columns)
    gene_list = list(mega_tbl.index.values)
    sample_info = pd.read_csv(args.table, sep="\t")
    for i in range(1, len(fname_list), 1):
        data = pd.read_csv(fname_list[i], sep="\t", index_col=0)
        mega_header.extend(list(data.columns))
        mega_tbl = pd.concat([mega_tbl, data])
    log_tf = np.log(np.array(mega_tbl))
    z_scored = stats.zscore(log_tf, axis = 1)
    mega_df = pd.DataFrame(z_scored, index=gene_list, columns=header)
    cbio_proj = sample_info.Cbio_project.unique()
    for proj in cbio_proj:
        new_fname = proj + ".rsem_merged_zscore.txt"
        cur = mega_df[[proj]]
        cur.fillna('NA', inplace=True)
        cur.to_csv(new_fname, sep='\t', mode='w', index=True)