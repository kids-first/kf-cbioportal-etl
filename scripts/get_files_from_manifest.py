#!/usr/bin/env python3

import argparse
import sys
import boto3
import concurrent.futures
import urllib3
import os
import pandas as pd
import sevenbridges as sbg
from sevenbridges.errors import SbgError
from sevenbridges.http.error_handlers import rate_limit_sleeper, maintenance_sleeper
import pdb


def download_aws(file_type):
    sub_file_list = list(selected.loc[selected["file_type"] == file_type, "s3_path"])
    for loc in sub_file_list:
        out = file_type + "/" + loc.split("/")[-1]
        parse_url = urllib3.util.parse_url(loc)
        try:
            dl_client.download_file(
                Bucket=parse_url.host, Key=parse_url.path.lstrip("/"), Filename=out
            )
        except Exception as e:
            sys.stderr.write(str(e) + " could not download from " + loc + "\n")
            sys.stderr.flush()
            err_types['aws download'] += 1
    sys.stderr.write("Completed downloading files for " + file_type + "\n")
    sys.stderr.flush()


def download_sbg(file_type):
    sub_file_list = list(selected.loc[selected["file_type"] == file_type, "file_id"])
    for loc in sub_file_list:
        try:
            sbg_file = api.files.get(loc)
        except Exception as e:
            sys.stderr.write('Failed to get file with id ' + loc + '\n')
            sys.stderr.write(str(e) + '\n')
            err_types['sbg get'] += 1
        out = file_type + "/" + sbg_file.name
        try:
            sbg_file.download(out)
        except Exception as e:
            err_types['sbg download'] += 1
            sys.stderr.write('Failed to download file with id ' + loc + '\n')
            sys.stderr.write(str(e) + '\n')
    return 0


def mt_type_download(file_type):
    """
    Download files from each desired file type at the same time
    """
    sys.stderr.write("Downloading " + file_type + " files\n")
    sys.stderr.flush()
    try:
        os.mkdir(file_type)
    except Exception as e:
        sys.stderr.write(
            str(e) + " error while making directory for " + file_type + "\n"
        )
    if args.profile:
        download_aws(file_type)
    else:
        download_sbg(file_type)


parser = argparse.ArgumentParser(description="Get all files for a project.")
parser.add_argument(
    "-m",
    "--manifest-list",
    action="store",
    dest="manifest",
    help="csv list of of genomic file location manifests",
)
parser.add_argument(
    "-f",
    "--file-types",
    action="store",
    dest="fts",
    help="csv list of workflow types to download",
)
parser.add_argument(
    "-p", "--profile", action="store", dest="profile", help="aws profile name. Leave blank if using sbg instead"
)
parser.add_argument(
    "-s", "--sbg-profile", action="store", dest="sbg_profile", help="sbg profile name. Leave blank if using AWS instead"
)
parser.add_argument(
    "-a", "--active-only", default=False, action="store_true", dest="active_only", help="Set to grab only active files. Recommended."
)


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
selected = manifest_concat[manifest_concat["file_type"].isin(file_types)]
# if active only flag set, further subset
if args.active_only:
    selected = selected[selected["status"] == 'active']
# remove vcfs as we only want mafs
pattern_del = ".vcf.gz"
filter = selected["file_name"].str.contains(pattern_del)
selected = selected[~filter]

out_file = "manifest_subset.tsv"
selected.to_csv(out_file, sep="\t", mode="w", index=False)

err_types = { 'aws download': 0, 'sbg get': 0, 'sbg download': 0 }
# download files by type
if args.profile is not None:
    session = boto3.Session(profile_name=args.profile)
    dl_client = session.client("s3")
elif args.sbg_profile is not None:
    config = sbg.Config(profile=args.sbg_profile)
    api = sbg.Api(config=config, error_handlers=[rate_limit_sleeper, maintenance_sleeper])
else:
    sys.stderr.write("Please set one of profile or sbg_profile\n")
    exit(1)

with concurrent.futures.ThreadPoolExecutor(16) as executor:
    results = {executor.submit(mt_type_download, ftype): ftype for ftype in file_types}

flag = 0
for key in err_types:
    if err_types[key]:
        flag += err_types[key]
        sys.stderr.write("ERROR: " + str(err_types[key]) + " " + key + " failure(s) occurred\n")
if flag:
    exit(1)
# for ftype in file_types:
#     mt_type_download(ftype)