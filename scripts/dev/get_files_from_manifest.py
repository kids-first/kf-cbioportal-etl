#!/usr/bin/env python3

import sys
import argparse
import concurrent.futures
import re
import pandas as pd
import pdb

parser = argparse.ArgumentParser(description='Get all files for a project.')
parser.add_argument('-m', '--manifest', action='store', dest='manifest', help='download manifest file csv list')
parser.add_argument('-w', '--wf-types', action='store', dest='wfs', help='csv list of workflow types to download')

args = parser.parse_args()
manifest = pd.read_csv(args.manifest)

wf_types = args.wfs.split(",")

selected = manifest[manifest['Workflow_type']].isin(wf_types)
# remove vcfs as we only want mafs
pattern_del = ".vcf.gz"
filter = selected['File_name'].str.contains(pattern_del)
selected = selected[~filter]

out_file = "manifest_subset.tsv"
selected.to_csv(expr_fname, sep='\t', mode='w', index=False)
