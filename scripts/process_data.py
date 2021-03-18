#!/usr/bin/env python3

import json
import sys
from pathlib import Path
import urllib.request
import os
import importlib
import requests
import base64
import pandas as pd

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
from .b_query_ds_by_bs_id import get_tumor_resources, query_dataservice_bs_id
from .c_create_data_sheets import create_master_dict
from .create_cbio_id_fname_tbl import process_ds_new

get_study_query = """query Study($id: ID!) {
  study(id: $id) {
    id
    name
    kfId
    description
    shortName
    projects {
      edges {
        node {
          name
          id
          projectId
        }
      }
    }
  }
}
"""


def get_diagnosis_cbio_disease_mapping(url):
    file = urllib.request.urlopen(url)
    d_dict = {}
    for line in file:
        (cbttc, cbio_short, cbio_long) = line.decode("utf-8").rstrip('\n').split('\t')
        d_dict[cbttc] = cbio_short
    return d_dict


def get_kf_id_to_cancer_type_mapping():
    d_dict = pd.read_csv("./REFS/kf_id_cancer_type_mapping.txt", delimiter='\t')
    return dict(d_dict.values)


def get_tumor_bs_mapped_normal_sample(config_data, resources):
    tumor_bs_mapped_normal = {}
    for resource in resources:
        if resource['atype'] == 'DNA':
            t_bs_id = resource.get('t_bs_id')
            n_bs_id = resource.get('n_bs_id')
            bs_data = query_dataservice_bs_id(n_bs_id, config_data['kf_url'], ['external_sample_id'], [], [], [])
            if bs_data['external_sample_id'] is None:
                sys.stderr.write('no normal data skip BS_ID: ' + n_bs_id + '\n')
            else:
                tumor_bs_mapped_normal[t_bs_id] = {
                    'specimen_id': n_bs_id,
                    'sample_id': format_smaple_id(config_data['dna_norm_id_style'], bs_data['external_sample_id'])
                }
    return tumor_bs_mapped_normal


def base64encode_study_id(study_id):
    encoded_bytes = base64.b64encode(("StudyNode:" + study_id).encode("utf-8"))
    return str(encoded_bytes, "utf-8")


if __name__ == '__main__':
    config_file = Path(__file__).parent / "../REFS/data_processing_config.json"
    config_data = {}
    with open(config_file) as f:
        config_data = json.load(f)

    if 'working_directory' not in config_data:
        parts = os.getcwd().split('/')
        config_data['working_directory'] = "/".join(parts[:-1])

    if config_data['working_directory'][-1] != '/':
        config_data['working_directory'] += '/'
    config_data['datasets'] = config_data['working_directory'] + 'datasets/'
    config_data['data_files'] = config_data['working_directory'] + 'data_files/'

    Path(config_data['datasets']).mkdir(parents=True, exist_ok=True)
    Path(config_data['data_files']).mkdir(parents=True, exist_ok=True)

    valid_extensions = []
    rna_ext_dict = config_data['rna_ext_list']
    for key in rna_ext_dict:
        valid_extensions.append(rna_ext_dict[key])
    dna_ext_dict = config_data['dna_ext_list']
    for key in dna_ext_dict:
        valid_extensions.append(dna_ext_dict[key])

    queried_study_ids = ['SD_M3DBXD12']
    encoded_study_ids = list(map(base64encode_study_id, queried_study_ids))

    kf_id_to_cancer_type_mapping = get_kf_id_to_cancer_type_mapping()

    for i in range(0, len(encoded_study_ids), 1):
        kf_study_id = queried_study_ids[i]
        encoded_study_id = encoded_study_ids[i]
        cancer_study_meta_info = {}
        cavatica_projects = []
        try:
            study_info = requests.post(
                config_data['kf_study_creator_url'],
                data=json.dumps({"query": get_study_query, "variables": {"id": encoded_study_id}}),
                headers={'Authorization': 'Bearer ' + config_data['kf_temp_bearer_token'],
                         'content_type': "application/json"}, )
            if study_info.json()['data'] is not None and study_info.json()['data']['study'] is not None:
                cancer_study_meta_info['cancer_study_identifier'] = kf_study_id
                cancer_study_meta_info['name'] = study_info.json()['data']['study']['name']
                cancer_study_meta_info['type_of_cancer'] = kf_id_to_cancer_type_mapping[kf_study_id]
                cancer_study_meta_info['description'] = study_info.json()['data']['study']['description']
                cancer_study_meta_info['short_name'] = study_info.json()['data']['study']['shortName']
                cancer_study_meta_info['groups'] = "PUBLIC"

                cavatica_project_id = 'kids-first-drc/' + (study_info.json()['data']['study']['kfId']).lower().replace(
                    "_", '-')
                study_projects = study_info.json()['data']['study']['projects']['edges']
                if len(study_projects) > 0:
                    for study_project in study_projects:
                        if cavatica_project_id == study_project['node']['projectId']:
                            cavatica_projects.append(study_project['node']['projectId'])
        except Exception as e:
            sys.stderr.write('Error ' + str(e) + ' occurred while trying to process ')
            sys.exit(1)

        filteredResources = get_resources_from_cavatica_projects(cavatica_projects, config_data)

        tum_out_fh = get_tumor_resources(filteredResources, config_data)
        tumor_bs_mapped_normal_sample = get_tumor_bs_mapped_normal_sample(filteredResources)

        dx_dict = get_diagnosis_cbio_disease_mapping(config_data['dx_tbl_fn'])
        create_master_dict(config_data, kf_study_id, tum_out_fh, dx_dict, tumor_bs_mapped_normal_sample)

        # download all files
        for filteredResource in filteredResources:
            for resource in filteredResource['resources']:
                resource.download(path=config_data['data_files']+resource.name)
        process_ds_new(config_data, cancer_study_meta_info, filteredResources)
