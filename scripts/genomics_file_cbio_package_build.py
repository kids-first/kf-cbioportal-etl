"""
This is a genomics file + clinical data file to cBio package conversion workflow script.
It is designed to work with standard KF somatic workflow outputs as well as DGD outputs.
Clinical files should have been produced ahead of time, while supporting sample ID file manifest, case_meta config json file, and data json config.
See README for prequisite details.
"""
#!/usr/bin/env python3
import sys
import os
import argparse
import json
import subprocess

parser = argparse.ArgumentParser(description='Download files (if needed), collate geomic files, organize load package.')
parser.add_argument('-t', '--table', action='store', dest='table',
                    help='Table with cbio project, kf bs ids, cbio IDs, and file names')
parser.add_argument('-m', '--manifest', action='store', dest='manifest', help='Download file manifest, if needed')
parser.add_argument('-c', '--cbio-config', action='store', dest='cbio_config', help='cbio case and meta config file')
parser.add_argument('-d', '--data-config', action='store', dest='data_config_file', help='json config file with data types and '
                                                                               'data locations')
parser.add_argument('-f', '--dgd-status', action='store', dest='dgd_status', help='Flag to determine load will have pbta/kf + dgd(both), kf/pbta only(kf), dgd-only(dgd)',
default='both', const='both', nargs='?', choices=['both', 'kf', 'dgd'])


args = parser.parse_args()
with open(args.data_config_file) as f:
    config_data = json.load(f)


def process_kf_maf():
    """
    Collate and process pbta/kf style mafs
    """
    status = 0
    return status


def process_append_dgd_maf():
    """
    Append DGD mafs to existing collated kf maf
    """
    status = 0
    return status


def process_dgd_maf():
    """
    Collate and process dgd style mafs
    """
    status = 0
    return status


def process_cnv():
    """
    Add gene info to CNV calls, merge into table, and create GISTIC-style output
    """
    status = 0
    return status


def process_seg():
    """
    Collate and process CNV seg files
    """
    status = 0
    return status


def process_rsem():
    """
    Merge rsem results by FPKM, calculate z-scores
    """
    status = 0
    return status


def process_kf_fusion():
    """
    Collate and process annoFuse output
    """
    status = 0
    return status


def append_dgd_fusion():
    """
    Append dgd fusions fo kf collated output
    """
    status = 0
    return status


def process_dgd_fusion():
    """
    Collated process DGD fusion output
    """
# This part is commented out for now as there is no good solution yet for files stored in different account - need to be able to run chopaws to get key
# if(args.manifest != None):
#     sys.stderr.write('Download manifest given. Downloading files\n')
#     dl_file_cmd = config_data['script_dir'] + 'get_files_from_manifest.py -m ' + args.manifest + ' -f ' + join(config_data['dl_file_type_list'] + ' -p saml')
#     subprocess.call(dl_file_cmd, shell=True)

with open(args.cbio_config) as f:
    config_meta_case = json.load(f)
for key in config_meta_case:
    if key.starts_with('meta_'):
        