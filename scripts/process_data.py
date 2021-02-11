#!/usr/bin/env python3

import sevenbridges as sbg
from sevenbridges.models.project import Project
from sevenbridges.models.file import File
from typing import List
import json
import sys
from requests import request
import argparse
import concurrent.futures
from time import sleep
import pdb
from pathlib import Path
import urllib.request
import os
import math
import re
import importlib

if __name__ == "__main__" and (__package__ is None or __package__ == ''):
    # replace the script's location in the Python search path by the main
    # scripts/ folder, above it, so that the importer package folder is in
    # scope and *not* directly in sys.path; see PEP 395
    sys.path[0] = str(Path(sys.path[0]).resolve().parent)
    __package__ = 'scripts'
    # explicitly load the package, which is needed on CPython 3.4 because it
    # doesn't include https://github.com/python/cpython/pull/2639
    importlib.import_module(__package__)

from .sample_id_builder_helper import format_smaple_id
from .a_merge_cavatica_manifests import get_resources_from_cavatica_projects
from .b_query_ds_by_bs_id import get_resource_information
from .c_create_data_sheets import create_master_dict
from .create_cbio_id_fname_tbl import process_ds_new

# START: STEP 1

def get_diagnosis_cbio_disease_mapping(url):
    file = urllib.request.urlopen(url)
    d_dict= {}
    for line in file:
        (cbttc, cbio_short, cbio_long) = line.decode("utf-8").rstrip('\n').split('\t')
        d_dict[cbttc] = cbio_short
    return d_dict

def get_norm_info(normal_objects):
    dna_samp_def = config_data['dna_norm_id_style']
    n_dict = {}
    for normal_object in normal_objects:
        if normal_object.get('BS_ID') is None or normal_object.get('external_sample_id') is None:
            sys.stderr.write('no normal data skip BS_ID: '+normal_object.get('BS_ID')+ '\n')
        else:
            n_dict[normal_object['BS_ID']] = format_smaple_id(dna_samp_def, normal_object['external_sample_id'])
    return n_dict

def get_tn_pair_data(resources):
    dna_pairs = {}
    for resource in resources:
        if resource['atype'] == 'DNA':
            t_bs_id = resource.get('t_bs_id')
            n_bs_id = resource.get('n_bs_id')
            if t_bs_id not in dna_pairs:
                dna_pairs[t_bs_id] = n_bs_id
            else:
                sys.stderr.write('WARN: tumor bs id ' + t_bs_id + ' already associated with a normal sample.  Skipping!\n')
                bs_ids_blacklist[t_bs_id] = True
    return dna_pairs

config_data =   {}

if __name__ == '__main__':
    config_file = Path(__file__).parent / "../REFS/data_processing_config.json"
    with open(config_file) as f:
        config_data = json.load(f)

    if 'working_directory' not in config_data:
        parts = os.getcwd().split('/')
        config_data['working_directory'] = "/".join(parts[:-1])
    
    if config_data['working_directory'][-1] != '/':
        config_data['working_directory'] += '/'
    config_data['datasets'] = config_data['working_directory'] + 'datasets/'
    config_data['data_files'] = config_data['working_directory'] + 'data_files/'

    valid_extensions = []
    rna_ext_dict = config_data['rna_ext_list']
    for key in rna_ext_dict:
        valid_extensions.append(rna_ext_dict[key])
    dna_ext_dict = config_data['dna_ext_list']
    for key in dna_ext_dict:
        valid_extensions.append(dna_ext_dict[key])

    pt_head = '\n'.join(config_data['ds_pt_desc'])

    # IMPORTANT! will use external sample id as sample id, and bs id as a specimen id
    samp_head = '\n'.join(config_data['ds_samp_desc'])

    filteredResources = get_resources_from_cavatica_projects(['zhangb1/pnoc003-cbio-data'], config_data, {})
    dna_pairs = get_tn_pair_data(filteredResources)

    resourceSet = {}
    # tum_out_fh, norm_out_fh = get_resource_information(filteredResources, config_data)

    # dx_dict = get_diagnosis_cbio_disease_mapping(config_data['dx_tbl_fn'])
    # norm_samp_id = get_norm_info(norm_out_fh)
    # create_master_dict(config_data, tum_out_fh, dx_dict, norm_samp_id, dna_pairs)
        
    #print(filteredResources)

    filteredResourceSet = {}

    for filteredResource in filteredResources:
        if filteredResource['atype'] == 'DNA':
            key = filteredResource['t_bs_id'] + filteredResource['n_bs_id'];
            filteredResourceSet[key] = filteredResource
        elif filteredResource['atype'] == 'RNA':
            filteredResourceSet[filteredResource['t_bs_id']] = filteredResource

    #download all files
    # out_dir = '/Users/kalletlak/Documents/temorary/data'

    # try:
    #     os.makedirs(out_dir)
    # except:
    #     sys.stderr.write(out_dir + ' already exists.\n')
    # for filteredResource in filteredResources:
    #     resource:File
    #     for resource in filteredResource['resources']:
    #         resource.download(path='/Users/kalletlak/Documents/temorary/data/'+resource.name)

    process_ds_new(config_data, filteredResourceSet)



