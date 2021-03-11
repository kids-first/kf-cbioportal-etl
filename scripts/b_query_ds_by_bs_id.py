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

bs_attrs = ['external_aliquot_id', 'analyte_type', 'source_text_tissue_type', 'source_text_tumor_descriptor',
            'composition', 'external_sample_id']
pt_attrs = ['external_id', 'gender', 'ethnicity', 'race']
outcome_attrs = ['age_at_event_days', 'vital_status']
dx_attrs = ['source_text_diagnosis', 'source_text_tumor_location', 'age_at_event_days']

def mt_query_calls(resource, url, norm_seen):
    try:
        tum_bs_id = resource.get('t_bs_id')
        norm_bs_id = resource.get('n_bs_id')
        norm_object = {}
        tumor_object = {}
        if norm_bs_id is not None and norm_bs_id not in norm_seen:
            norm_object = query_dataservice_bs_id(norm_bs_id, bs_attrs, pt_attrs, outcome_attrs, url)
            if bool(norm_object):
                norm_object['BS_ID'] = norm_bs_id
            norm_seen[norm_bs_id] = 1
        tumor_object = query_dataservice_bs_id(tum_bs_id, bs_attrs, pt_attrs, outcome_attrs, url)
        if bool(tumor_object):
            tumor_object['BS_ID'] = tum_bs_id
        return tumor_object, norm_object
    except Exception as e:
        sys.stderr.write('Error ' + str(e) + ' occurred while trying to process ')
        exit(1)


def query_dataservice_bs_id(bs_id, bs_attrs, pt_attrs, outcome_attrs, url):
    result = {}
    try:
        bs_url = url + '/biospecimens/' + bs_id
        bs_info = request('GET', bs_url)
        if bs_info.json()['_status']['code'] == 404:
            sys.stderr.write(bs_id + ' not found!\n')
            sys.stderr.flush()
            return result
        if not bs_info.json()['results']['visible']:
            sys.stderr.write(bs_id + ' was set to hide.  skipping!\n')
            sys.stderr.flush()
            return result

        for attr in bs_attrs:
            result[attr] = bs_info.json()['results'][attr]

         # fetch diagnosis
        try:
            result['diagnosis'] = []
            dx_url = url + bs_info.json()['_links']['diagnoses']
            dx_obj = request('GET', dx_url)
            for cur_res in dx_obj.json()['results']:
                values = {}
                for attr in dx_attrs:
                    values[attr] = cur_res[attr]
                result['diagnosis'].append(values)

        except Exception as e:
            sys.stderr.write('Got error ' + str(e) + ',could not retrieve info from ' + dx_url + '\n')

        # fetch patient information
        try:
            pt_url = bs_info.json()['_links']['participant']
            pt_info = request('GET', url + pt_url)
            result['PT_ID'] = pt_info.json()['results']['kf_id']
            for attr in pt_attrs:
                result[attr] = pt_info.json()['results'][attr]
        except Exception as e:
            sys.stderr.write('Got error ' + str(e) + ', could not retrieve info from ' + pt_url + '\n')

        # fetch patient outcomes
        outcome_url = pt_info.json()['_links']['outcomes']
        try:
            outcome_info = request('GET', url + outcome_url)
            # default to last outcome in list
            try:
                if len(outcome_info.json()['results']):
                    for attr in outcome_attrs:
                        # had to add outcome_ prefix, there are duplicate attr from dx_attrs. expecially age_at_event_days
                        result['outcome_'+attr] = outcome_info.json()['results'][-1][attr]
                else:
                    for attr in outcome_attrs:
                        result['outcome_'+attr] = None
                        sys.stderr.write("WARN: " + pt_info.json()['results']['kf_id'] + " has no outcome data\n")
            except Exception as e:
                sys.stderr.write(str(e) + "\n")

        except Exception as e:
            sys.stderr.write('Got error ' + str(e) + ', could not retrieve info from ' + outcome_url + '\n')

    except Exception as e:
        sys.stderr.write('Got error ' + str(e) + ',could not retrieve info from ' + bs_url + '\n')

    return result

def get_resource_information(resources, config_data):
    norm_out_fh = []
    tum_out_fh = []
    norm_seen = {}
    with concurrent.futures.ThreadPoolExecutor(config_data['threads']) as executor:
        results = {executor.submit(mt_query_calls, resource, config_data['kf_url'], norm_seen): resource for resource in resources}
        for result in concurrent.futures.as_completed(results):
            (tum_info, norm_info) = result.result()
            if bool(tum_info):
                tum_out_fh.append(tum_info)
            if bool(norm_info):
                norm_out_fh.append(norm_info)
    return tum_out_fh,norm_out_fh
