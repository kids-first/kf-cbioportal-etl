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

def get_filtered_resources(resources:List[File], config_data, skip_dict = {}):

    rna_ext_dict = config_data['rna_ext_list']
    rna_ext_list = []
    for key in rna_ext_dict:
        rna_ext_list.append(rna_ext_dict[key])
    dna_ext_dict = config_data['dna_ext_list']
    dna_ext_list = []
    for key in dna_ext_dict:
        dna_ext_list.append(dna_ext_dict[key])

    cav_dict = {}
    for resource in resources:
        fname = resource.name
        t_bs_id = resource.metadata['Kids First Biospecimen ID Tumor']
        if t_bs_id is None:
            t_bs_id = resource.metadata['Kids First Biospecimen ID']
        project = resource.project
        n_bs_id = resource.metadata['Kids First Biospecimen ID Normal']

        if t_bs_id in skip_dict or n_bs_id in skip_dict:
            sys.stderr.write('BS ID in blacklist.  Skipping ' + t_bs_id)
            continue

        parts = fname.split('.')
        ext = ".".join(parts[1:])

        if ext in rna_ext_list:
            if t_bs_id not in cav_dict:
                cav_dict[t_bs_id] = {}
                cav_dict[t_bs_id]['atype'] = "RNA"
                cav_dict[t_bs_id]['fname'] = []
                cav_dict[t_bs_id]['resources'] = []
                cav_dict[t_bs_id]['t_bs_id'] = t_bs_id
                cav_dict[t_bs_id]['project'] = []
            cav_dict[t_bs_id]['project'].append(project)
            cav_dict[t_bs_id]['fname'].append(fname)
            cav_dict[t_bs_id]['resources'].append(resource)
            
        elif ext in dna_ext_list:
            dkey = t_bs_id + n_bs_id
            if dkey not in cav_dict:
                cav_dict[dkey] = {}
                cav_dict[dkey]['atype'] = "DNA"
                cav_dict[dkey]['fname'] = []
                cav_dict[dkey]['resources'] = []
                cav_dict[dkey]['t_bs_id'] = t_bs_id
                cav_dict[dkey]['n_bs_id'] = n_bs_id
                cav_dict[dkey]['project'] = []
            cav_dict[dkey]['project'].append(project)
            cav_dict[dkey]['fname'].append(fname)
            cav_dict[dkey]['resources'].append(resource)
    bs_ids_blacklist = {}
    dna_pairs = {}
    for key, value in cav_dict.items():
        if value['atype'] == 'DNA':
            t_bs_id = value.get('t_bs_id')
            n_bs_id = value.get('n_bs_id')
            if t_bs_id not in dna_pairs:
                dna_pairs[t_bs_id] = n_bs_id
            else:
                sys.stderr.write('WARN: tumor bs id ' + t_bs_id + ' already associated with a normal sample.  Skipping!\n')
                bs_ids_blacklist[t_bs_id] = True
    filteredResources  = []
    for key, value in cav_dict.items():
        t_bs_id = value['t_bs_id']
        if t_bs_id is not None and t_bs_id not in bs_ids_blacklist:
            filteredResources.append(value)
    return filteredResources


def get_resources_from_cavatica_projects(project_ids, config_data, skip_dict = {}):
    valid_extensions = []
    rna_ext_dict = config_data['rna_ext_list']
    for key in rna_ext_dict:
        valid_extensions.append(rna_ext_dict[key])
    dna_ext_dict = config_data['dna_ext_list']
    for key in dna_ext_dict:
        valid_extensions.append(dna_ext_dict[key])

    api = sbg.Api(url='https://cavatica-api.sbgenomics.com/v2', token='296f647e655b4ff2acbc06f92a56b733')
    resources: List[File] =[]
    for project_id in project_ids:
        project:Project = api.projects.get(id=project_id)
        #TODO: multiple project
        #projects = api.projects.query(owner='zhangb1', name='PNOC003 Cbio Data')
        folders: List[File] = list(project.get_files())
        for folder in folders:
            if folder.is_folder:
                for x in range(0, 5, 1):
                    response = api.files.bulk_get(folder.list_files(offset=x, limit=1))
                    for record in response:
                        if record.valid:
                            name = record.resource.name
                            parts = name.split('.')
                            ext = ".".join(parts[1:])
                            if ext in valid_extensions:
                                
                                resources.append(record.resource)
                        else:
                            print(record.error)
    return get_filtered_resources(resources, config_data, skip_dict)