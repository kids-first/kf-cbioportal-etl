#!/usr/bin/env python3

import argparse
import pdb
import sys
import boto3
import concurrent.futures
import urllib3
import os
import pandas as pd


def mt_type_download(file_type):
    """
    Download files from each desired file type at the same time
    """
    sys.stderr.write("Downloading " + file_type + " files\n")
    sys.stderr.flush()
    sub_file_list = list(selected.loc[selected['file_type'] == file_type, 's3_path'])
    try:
        os.mkdir(file_type)
    except Exception as e:
        sys.stderr.write(str(e) + " error while making directory for " + file_type + "\n")
    for loc in sub_file_list:
        out = file_type + "/" + loc.split("/")[-1]
        parse_url = urllib3.util.parse_url(loc)
        try:
            dl_client.download_file(Bucket=parse_url.host, Key=parse_url.path.lstrip("/"), Filename=out)
        except Exception as e:
            sys.stderr.write(str(e) + " could not download from " + loc + "\n")
            sys.stderr.flush()
    sys.stderr.write("Completed downloading files for " + file_type + "\n")
    sys.stderr.flush()


parser = argparse.ArgumentParser(description='Get all files for a project.')
parser.add_argument('-m', '--manifest-list', action='store', dest='manifest', help='csv list of of genomic file location manifests')
parser.add_argument('-f', '--file-types', action='store', dest='fts', help='csv list of workflow types to download')
parser.add_argument('-p', '--profile', action='store', dest='profile', help='aws profile name')

args = parser.parse_args()
# concat multiple possible manifests
sys.stderr.write("Concatenating manifests\n")
sys.stderr.flush()
manifest_list = args.manifest.split(",")
manifest_concat = pd.DataFrame()
for manifest in manifest_list:
    sys.stderr.write("Processing " + manifest + "\n")
    current = pd.read_csv(manifest, sep=None)
    if manifest_concat.empty:
        manifest_concat = current.copy()
    else:
        manifest_concat = manifest_concat.append(current, ignore_index=True)

file_types = args.fts.split(",")
# subset concatenated manifests
sys.stderr.write("Subsetting concatenated manifest\n")
sys.stderr.flush()
# change col names to lower case for input compatibility
manifest_concat.columns = manifest_concat.columns.str.lower()
selected = manifest_concat[manifest_concat['file_type'].isin(file_types)]
# remove vcfs as we only want mafs
pattern_del = ".vcf.gz"
filter = selected['file_name'].str.contains(pattern_del)
selected = selected[~filter]

out_file = "manifest_subset.tsv"
selected.to_csv(out_file, sep='\t', mode='w', index=False)

# download files by type
session = boto3.Session(profile_name=args.profile)
dl_client = session.client('s3')
# pdb.set_trace()

with concurrent.futures.ThreadPoolExecutor(16) as executor:
    results = {executor.submit(mt_type_download, ftype): ftype for ftype in file_types}

