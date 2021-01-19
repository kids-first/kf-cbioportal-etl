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
            norm_object = query_dataservice_bs_id(norm_bs_id, bs_attrs, pt_attrs, dx_attrs, outcome_attrs, url)
            if bool(norm_object):
                norm_object['BS_ID'] = norm_bs_id
            norm_seen[norm_bs_id] = 1
        tumor_object = query_dataservice_bs_id(tum_bs_id, bs_attrs, pt_attrs, dx_attrs, outcome_attrs, url)
        if bool(tumor_object):
            tumor_object['BS_ID'] = tum_bs_id
        return tumor_object, norm_object
    except Exception as e:
        sys.stderr.write('Error ' + str(e) + ' occurred while trying to process ')
        exit(1)


def query_dataservice_bs_id(bs_id, bs_attrs, pt_attrs, dx_attrs, outcome_attrs, url):
    bs_url = url + '/biospecimens/' + bs_id
    # sys.stderr.write(bs_url + '\n')
    bs_info = None
    try:
        bs_info = request('GET', bs_url)
    except Exception as e:
        sys.stderr.write('Got error ' + str(e) + ',could not retrieve info from ' + bs_url + ', pausing and retrying\n')
        sys.stderr.flush()
        sleep(1)
        bs_info = request('GET', bs_url)
    result = {}
    if bs_info.json()['_status']['code'] == 404:
        # result.append(bs_info.json()['_status']['message'])
        sys.stderr.write(bs_id + ' not found!\n')
        sys.stderr.flush()
        return result
    elif not bs_info.json()['results']['visible']:
        sys.stderr.write(bs_id + ' was set to hide.  skipping!\n')
        sys.stderr.flush()
        return result
    dx_url = url + bs_info.json()['_links']['diagnoses']
    dx_dict = {}
    dx_obj = 'NoDX'
    try:
        dx_obj = request('GET', dx_url) if len(dx_attrs) > 0 else 'NoDX'
    except Exception as e:
        sys.stderr.write('Got error ' + str(e) + ',could not retrieve info from ' + dx_url + ', pausing and retrying\n')
        sys.stderr.flush()
        sleep(1)
        dx_obj = request('GET', dx_url) if len(dx_attrs) > 0 else 'NoDX'
    # dir(bs_info)
    pt_url = bs_info.json()['_links']['participant']
    pt_info = None
    try:
        pt_info = request('GET', url + pt_url)
    except Exception as e:
        sys.stderr.write('Got error ' + str(e) + ', could not retrieve info from ' + pt_url + ', pausing and retrying\n')
        sys.stderr.flush()
        sleep(1)
        pt_info = request('GET', url + pt_url)
    result['PT_ID'] = pt_info.json()['results']['kf_id']
    for attr in bs_attrs:
        result[attr] = bs_info.json()['results'][attr]
    for attr in pt_attrs:
        result[attr] = pt_info.json()['results'][attr]
    for attr in dx_attrs:
        values = []
        for cur_res in dx_obj.json()['results']:
            for dx_split in cur_res['source_text_diagnosis'].split(','):
                if attr == 'source_text_diagnosis':
                    values.append(dx_split.strip()) 
                else:   
                    values.append(str(cur_res[attr]))
        result[attr] = ';'.join(values)

    outcome_url = pt_info.json()['_links']['outcomes']
    outcome_info = None
    try:
        outcome_info = request('GET', url + outcome_url)
    except Exception as e:
        sys.stderr.write('Got error ' + str(e) + ', could not retrieve info from ' + outcome_url + ', pausing and retrying\n')
        sys.stderr.flush()
        sleep(1)
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
    return result

def get_resource_information(resources, config_data):
    x = 1
    m = 100
    norm_out_fh = []
    tum_out_fh = []
    norm_seen = {}
    with concurrent.futures.ThreadPoolExecutor(config_data['threads']) as executor:

        results = {executor.submit(mt_query_calls, resource, config_data['kf_url'], norm_seen): resource for resource in resources}
        for result in concurrent.futures.as_completed(results):
            if x % m == 0:
                sys.stderr.write('Processed ' + str(x) + ' lines\n')
                sys.stderr.flush()
            (tum_info, norm_info) = result.result()
            if bool(tum_info):
                tum_out_fh.append(tum_info)
            if bool(norm_info):
                norm_out_fh.append(norm_info)
            x += 1
    return tum_out_fh,norm_out_fh
