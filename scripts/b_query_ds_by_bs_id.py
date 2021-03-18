#!/usr/bin/env python3

import sys
from requests import request
import concurrent.futures

bs_attrs = ['external_aliquot_id', 'analyte_type', 'source_text_tissue_type', 'source_text_tumor_descriptor',
            'composition', 'external_sample_id']
pt_attrs = ['external_id', 'gender', 'ethnicity', 'race']
outcome_attrs = ['age_at_event_days', 'vital_status']
dx_attrs = ['source_text_diagnosis', 'source_text_tumor_location', 'age_at_event_days']


def query_dataservice_bs_id(bs_id, url, bs_attrs, dx_attrs, pt_attrs, outcome_attrs):
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
        if len(dx_attrs) > 0:
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
                        result['outcome_' + attr] = outcome_info.json()['results'][-1][attr]
                else:
                    for attr in outcome_attrs:
                        result['outcome_' + attr] = None
                        sys.stderr.write("WARN: " + pt_info.json()['results']['kf_id'] + " has no outcome data\n")
            except Exception as e:
                sys.stderr.write(str(e) + "\n")

        except Exception as e:
            sys.stderr.write('Got error ' + str(e) + ', could not retrieve info from ' + outcome_url + '\n')

    except Exception as e:
        sys.stderr.write('Got error ' + str(e) + ',could not retrieve info from ' + bs_url + '\n')
    result['BS_ID'] = bs_id
    return result


def get_tumor_resources(resources, config_data):
    tum_out_fh = []
    with concurrent.futures.ThreadPoolExecutor(config_data['threads']) as executor:
        results = {
            executor.submit(query_dataservice_bs_id, resource.get('t_bs_id'), config_data['kf_url'], bs_attrs, dx_attrs,
                            pt_attrs, outcome_attrs): resource for resource in resources}
        for result in concurrent.futures.as_completed(results):
            tum_info = result.result()
            if bool(tum_info):
                tum_out_fh.append(tum_info)
    return tum_out_fh
