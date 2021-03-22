#!/usr/bin/env python3

import json
import sys
from pathlib import Path
import urllib.request
import requests
import base64
import pandas as pd
import argparse

from cavatica_resources import get_resources_from_cavatica_projects
from cbioportal_resources import process_study_resources

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


def base64encode_study_id(study_id):
    encoded_bytes = base64.b64encode(("StudyNode:" + study_id).encode("utf-8"))
    return str(encoded_bytes, "utf-8")


def main():
    parser = argparse.ArgumentParser()
    # parser.add_argument("--kf_ids", help='List of kf study ids')
    parser.add_argument("-i", "--kf_ids", nargs='+', type=str, required=True, help='List of kf study ids')
    parser.add_argument("-c", "--cavatica_token", type=str, required=True, help='Cavatica token')
    parser.add_argument("-k", "--kf_token", type=str, required=True, help='KF bearer token token')
    args = parser.parse_args()

    config_file = Path(__file__).parent / "../REFS/data_processing_config.json"
    config_data = {}
    with open(config_file) as f:
        config_data = json.load(f)

    config_data['kf_bearer_token'] = args.kf_token
    config_data['cavatica_token'] = args.cavatica_token

    config_data['datasets'] = str(Path.home()) + '/cbioportal/datasets/'
    config_data['data_files'] = str(Path.home()) + '/cbioportal/data_files/'

    Path(config_data['datasets']).mkdir(parents=True, exist_ok=True)
    Path(config_data['data_files']).mkdir(parents=True, exist_ok=True)

    valid_extensions = []
    rna_ext_dict = config_data['rna_ext_list']
    for key in rna_ext_dict:
        valid_extensions.append(rna_ext_dict[key])
    dna_ext_dict = config_data['dna_ext_list']
    for key in dna_ext_dict:
        valid_extensions.append(dna_ext_dict[key])

    encoded_study_ids = list(map(base64encode_study_id, args.kf_ids))

    kf_id_to_cancer_type_mapping = get_kf_id_to_cancer_type_mapping()
    dx_dict = get_diagnosis_cbio_disease_mapping(config_data['dx_tbl_fn'])

    for i in range(0, len(encoded_study_ids), 1):
        kf_study_id = args.kf_ids[i]
        encoded_study_id = encoded_study_ids[i]
        cancer_study_meta_info = {}
        cavatica_projects = []
        try:
            study_info = requests.post(
                config_data['kf_study_creator_url'],
                data=json.dumps({"query": get_study_query, "variables": {"id": encoded_study_id}}),
                headers={'Authorization': 'Bearer ' + config_data['kf_bearer_token'],
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
                            break
        except Exception as e:
            sys.stderr.write('Error ' + str(e) + ' occurred while trying to process\n')
            sys.exit(1)

        resources = get_resources_from_cavatica_projects(cavatica_projects, config_data)
        process_study_resources(config_data, cancer_study_meta_info, resources, dx_dict)


if __name__ == '__main__':
    main()
