#!/usr/bin/env python3

import argparse
import pdb
import sys

import pandas as pd

parser = argparse.ArgumentParser(description='Get all files for a project.')
parser.add_argument('-m', '--manifest-list', action='store', dest='manifest', help='csv list of of genomic file location manifests')
parser.add_argument('-f', '--file-types', action='store', dest='fts', help='csv list of workflow types to download')

args = parser.parse_args()
# concat multiple possible manifests
sys.stderr.write("Concatenating manifests\n")
manifest_list = args.manifest.split(",")
manifest_concat = pd.DataFrame()
for manifest in manifest_list:
    sys.stderr.write("Processing " + manifest + "\n")
    current = pd.read_csv(manifest)
    if manifest_concat.empty:
        manifest_concat = current.copy()
    else:
        manifest_concat = manifest_concat.append(current, ignore_index=True)

file_types = args.fts.split(",")
# subset concatenated manifests
sys.stderr.write("Subsetting concatenated manifest\n")
selected = manifest_concat[manifest_concat['File_type'].isin(file_types)]
# remove vcfs as we only want mafs
pattern_del = ".vcf.gz"
filter = selected['File_name'].str.contains(pattern_del)
selected = selected[~filter]

out_file = "manifest_subset.tsv"
selected.to_csv(out_file, sep='\t', mode='w', index=False)


