#!/usr/bin/env python3

import argparse
import sys
import boto3
import botocore
import concurrent.futures
import urllib3
import os
import pandas as pd
import sevenbridges as sbg
from sevenbridges.http.error_handlers import rate_limit_sleeper, maintenance_sleeper
from time import sleep
import pdb


def download_aws(file_type):
    """
    Function to download from AWS to file_type dir if AWS table was given
    """
    for key in key_dict:
        current = key_dict[key]["manifest"]
        dl_client = key_dict[key]["dl_client"]
        sub_file_list = list(current.loc[current["file_type"] == file_type, "s3_path"])
        print("Grabbing {} files using key {}".format(len(sub_file_list), key), file=sys.stderr)
        for loc in sub_file_list:
            out = file_type + "/" + loc.split("/")[-1]
            if args.overwrite or not os.path.isfile(out):
                parse_url = urllib3.util.parse_url(loc)
                try:
                    dl_client.download_file(
                        Bucket=parse_url.host, Key=parse_url.path.lstrip("/"), Filename=out
                    )
                except Exception as e:
                    print("{} could not download from {} using key {}".format(e, loc, key), file=sys.stderr)
                    sys.stderr.flush()
                    err_types['aws download'] += 1
            else:
                sys.stderr.write("Skipping " + out + " it exists and overwrite not set\n")
    print("Completed downloading files for " + file_type, file=sys.stderr)
    sys.stderr.flush()


def download_sbg(file_type):
    """
    Function to download from SBG to file_type dir if SBG profile was given
    """
    sub_file_list = list(selected.loc[selected["file_type"] == file_type, "file_id"])
    # Sort of a "trust fall" that if aws bucket exists, skip SBG download
    if args.aws_tbl:
        sub_file_list = list(selected.loc[(selected["file_type"] == file_type) & (selected["s3_path"].isna()), "file_id"])
    for loc in sub_file_list:
        try:
            sbg_file = api.files.get(loc)
        except Exception as e:
            print("Failed to get file with id " + loc + ". Will try once more in 3 seconds", file=sys.stderr)
            print(e, file=sys.stderr)
            try:
                sleep(3)
                sbg_file = api.files.get(loc)
                print("Success on second try!", file=sys.stderr)
            except Exception as e:
                print("Failed on second attempt", file=sys.stderr)
                err_types['sbg get'] += 1
        out = file_type + "/" + sbg_file.name
        if args.overwrite or not os.path.isfile(out):
            try:
                sbg_file.download(out)
            except Exception as e:
                err_types['sbg download'] += 1
                print("Failed to download file with id " + loc, file=sys.stderr)
                print(e, file=sys.stderr)
        else:
            print("Skipping " + out + " it exists and overwrite not set", file=sys.stderr)
    return 0


def mt_type_download(file_type):
    """
    Download files from each desired file type at the same time
    """
    if len(selected.loc[selected["file_type"] == file_type]):
        print("Downloading " + file_type + " files", file=sys.stderr)
        try:
            os.mkdir(file_type)
        except Exception as e:
            print("{} error while making directory for {}".format(e, file_type), file=sys.stderr)
        if args.aws_tbl:
            download_aws(file_type)
        if args.sbg_profile:
            download_sbg(file_type)
    else:
        print("WARN: No files of type " + file_type + " in which file_id and s3_path is not NA. Skipping!", file=sys.stderr)
    sys.stderr.flush()


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
    "-t", "--aws-tbl", action="store", dest="aws_tbl", help="Table with bucket name and keys to subset on"
)
parser.add_argument(
    "-s", "--sbg-profile", action="store", dest="sbg_profile", help="sbg profile name. Leave blank if using AWS instead"
)
parser.add_argument(
    "-c", "--cbio", action="store", dest="cbio", help="Add cbio manifest to limit downloads"
)
parser.add_argument(
    "-a", "--active-only", default=False, action="store_true", dest="active_only", help="Set to grab only active files. Recommended."
)
parser.add_argument(
    "-r", "--rm-na", default=False, action="store_true", dest="rm_na", help="Remove entries where file_id and s3_path are NA. \
    Useful for studies (like pbta_all) with external files not to be downloaded when using cbio_file_name input file as manifest"
)
parser.add_argument(
    "-d", "--debug", default=False, action="store_true", dest="debug", help="Just output manifest subset to see what would be grabbed"
)
parser.add_argument(
    "-o", "--overwrite", default=False, action="store_true", dest="overwrite", help="If set, overwrite if file exists"
)


args = parser.parse_args()
# concat multiple possible manifests
print("Concatenating manifests", file=sys.stderr)
sys.stderr.flush()
manifest_list = args.manifest.split(",")
manifest_df_list = []
for manifest in manifest_list:
    sys.stderr.write("Processing " + manifest + "\n")
    manifest_df_list.append(pd.read_csv(manifest, sep=None, na_values=['NA']))
# In the event that s3_path is empty, replace with str to trigger later sbg download
manifest_concat = pd.concat(manifest_df_list, ignore_index=True)
# if using a cbio name file as manifest, drop conflicting column
if "File_Type" in manifest_concat.columns:
    del manifest_concat["File_Type"]
# change col names to lower case for input compatibility
manifest_concat.columns = manifest_concat.columns.str.lower()
# file_types is actually a requirement, so grab from table if not provided
if args.fts:
    file_types = args.fts.split(",")
else:
    print("No file types provided, using table values", file=sys.stderr)
    file_types = manifest_concat['file_type'].unique().tolist()
# subset concatenated manifests
print("Subsetting concatenated manifest", file=sys.stderr)
sys.stderr.flush()
selected = manifest_concat[manifest_concat["file_type"].isin(file_types)]
# if active only flag set, further subset
if args.active_only:
    print("active only flag given, will limit to that", file=sys.stderr)
    selected = selected[selected["status"] == 'active']
if args.rm_na:
    print("Drop na flag given. Will drop rows in which file_id and s3_path are NA", file=sys.stderr)
    selected = selected[ ~(selected['file_id'].isna() & selected['s3_path'].isna()) ]
# if cbio manifest given, limit to that
if args.cbio:
    print("cBio manifest provided, limiting downloads to matching IDs", file=sys.stderr)
    cbio_data = pd.read_csv(args.cbio, sep=None, na_values=['NA'])
    specimen_list = cbio_data.T_CL_BS_ID.unique()
    selected = selected[selected["biospecimen_id"].isin(specimen_list)]
# remove vcfs as we only want mafs
pattern_del = ".vcf.gz"
filter = selected["file_name"].str.contains(pattern_del)
selected = selected[~filter]

out_file = "manifest_subset.tsv"
selected.to_csv(out_file, sep="\t", mode="w", index=False)
if args.debug:
    sys.stderr.write('Debug flag given. No downloads actually happen, just a manifest subset to preview\n')
    exit(1)
err_types = { 'aws download': 0, 'sbg get': 0, 'sbg download': 0 }
# download files by type
check = 0
if args.aws_tbl is not None:
    check = 1
    # from https://stackoverflow.com/questions/66041582/connection-pool-is-full-warning-while-reading-s3-objects-via-multiple-threads
    client_config = botocore.config.Config(
    max_pool_connections=128
)
    # setting up a key dict that, for each aws key, has an associated sesion and manifest to download with
    key_dict = {}
    bucket_errs = 0
    with open(args.aws_tbl) as kl:
        for line in kl:
            (bucket, key) = line.rstrip('\n').split('\t')
            if key not in key_dict:
                key_dict[key] = {}
                key_dict[key]['manifest'] = selected[selected['s3_path'].str.contains(bucket)]
                key_dict[key]['session'] = boto3.Session(profile_name=key)
                key_dict[key]['dl_client'] = key_dict[key]['session'].client("s3", config=client_config)
            else:
                key_dict[key]['manifest'] = pd.concat([key_dict[key]['manifest'], selected[selected['s3_path'].str.startswith(bucket)]], ignore_index=True)
            # Test bucket access with that key, if it fails, print error then kill to not waste time
            parse_url = urllib3.util.parse_url(bucket)
            try:
                key_dict[key]['dl_client'].list_objects(Bucket=parse_url.host)
            except Exception as e:
                bucket_errs = 1
                print(e, file=sys.stderr)
                print("Bucket access ERROR: {}\t{}".format(bucket, key), file=sys.stderr)
    if bucket_errs:
        exit(1)
if args.sbg_profile is not None:
    check = 1
    config = sbg.Config(profile=args.sbg_profile)
    api = sbg.Api(config=config, error_handlers=[rate_limit_sleeper, maintenance_sleeper])
if not check:
    sys.stderr.write("Please provide at least one of aws_tbl or sbg_profile\n")
    exit(1)

with concurrent.futures.ThreadPoolExecutor(16) as executor:
    results = {executor.submit(mt_type_download, ftype): ftype for ftype in file_types}

# I leave this here for debugging
# for ftype in file_types:
#     mt_type_download(ftype)

flag = 0
for key in err_types:
    if err_types[key]:
        flag += err_types[key]
        print("ERROR: {} {} failure(s) occurred".format(err_types[key], key), file=sys.stderr)
if flag:
    exit(1)
else:
    print("Downloads complete. No errors caught", file=sys.stderr)
