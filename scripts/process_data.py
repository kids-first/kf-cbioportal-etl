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
# START: STEP 1

def get_diagnosis_cbio_disease_mapping(url):
    file = urllib.request.urlopen(url)
    d_dict= {}
    for line in file:
        (cbttc, cbio_short, cbio_long) = line.decode("utf-8").rstrip('\n').split('\t')
        d_dict[cbttc] = cbio_short
    return d_dict

def get_resources(resources:[File], config_data, skip_dict = {}):

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
                cav_dict[t_bs_id]['t_bs_id'] = t_bs_id
                cav_dict[t_bs_id]['project'] = []
            cav_dict[t_bs_id]['project'].append(project)
            cav_dict[t_bs_id]['fname'].append(fname)
        elif ext in dna_ext_list:
            dkey = t_bs_id + n_bs_id
            if dkey not in cav_dict:
                cav_dict[dkey] = {}
                cav_dict[dkey]['atype'] = "DNA"
                cav_dict[dkey]['fname'] = []
                cav_dict[dkey]['t_bs_id'] = t_bs_id
                cav_dict[dkey]['n_bs_id'] = n_bs_id
                cav_dict[dkey]['project'] = []
            cav_dict[dkey]['project'].append(project)
            cav_dict[dkey]['fname'].append(fname)
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
    return filteredResources,dna_pairs

def mt_query_calls(resource):
    try:
        tum_bs_id = resource.get('t_bs_id')
        norm_bs_id = resource.get('n_bs_id')
        norm_object = {}
        tumor_object = {}
        if norm_bs_id is not None and norm_bs_id not in norm_seen:
            norm_object = query_dataservice_bs_id(norm_bs_id, bs_attrs, pt_attrs, dx_attrs, outcome_attrs)
            if bool(norm_object):
                norm_object['BS_ID'] = norm_bs_id
            norm_seen[norm_bs_id] = 1
        tumor_object = query_dataservice_bs_id(tum_bs_id, bs_attrs, pt_attrs, dx_attrs, outcome_attrs)
        if bool(tumor_object):
            tumor_object['BS_ID'] = tum_bs_id

        return tumor_object, norm_object
    except Exception as e:
        sys.stderr.write('Error ' + str(e) + ' occurred while trying to process ')
        exit(1)

def get_norm_info(normal_objects):
    dna_samp_def = config_data['dna_norm_id_style']
    n_dict = {}
    for normal_object in normal_objects:
        if normal_object.get('BS_ID') is None or normal_object.get('external_sample_id') is None:
            sys.stderr.write('no normal data skip BS_ID: '+normal_object.get('BS_ID')+ '\n')
        else:
            n_dict[normal_object['BS_ID']] = format_smaple_id(dna_samp_def, normal_object['external_sample_id'])
    return n_dict

def query_dataservice_bs_id(bs_id, bs_attrs, pt_attrs, dx_attrs, outcome_attrs):
    url = config_data['kf_url']
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

def create_master_dict(tumor_objects, dx_dict, norm_samp_id):
    blacklist = {}
    master_dict = {}
    dna_samp_def = config_data['dna_samp_id_style']
    rna_samp_def = config_data['rna_samp_id_style']

    # collate tum info with DNA pair info and norm info by dx
    sys.stderr.write('Collating information for data sheet\n')
    for tumor_object in tumor_objects:
        bs_id = tumor_object['BS_ID']
        pt_id = tumor_object['PT_ID']
        ana_type = tumor_object['analyte_type']
        samp_id = format_smaple_id(dna_samp_def, tumor_object['external_sample_id'])
        samp_type = tumor_object['composition']
        # if samp_type == 'Derived Cell Line':
        #     samp_id += '-CL'
        #     if args.cl_supp is not None and bs_id in cl_supp:
        #         samp_id += '-' + cl_supp[bs_id]
        cdx_list = tumor_object['source_text_diagnosis'].split(';')
        loc_list = tumor_object['source_text_tumor_location'].split(';')
        age_list = []
        if tumor_object.get('age_at_event_days') is not None:
            age_list =  tumor_object.get('age_at_event_days').split(';')
        temp_dx = {}
        for i in range(0, len(cdx_list), 1):
            if cdx_list[i] == "N/A":
                cdx_list[i] = "Not Reported"
            if cdx_list[i] != '' and cdx_list[i] in dx_dict:
                cur_dx = dx_dict[cdx_list[i]]
                if cur_dx not in temp_dx:
                    temp_dx[cur_dx] = {}
                    temp_dx[cur_dx]['dx'] = []
                    temp_dx[cur_dx]['loc'] = []
                    temp_dx[cur_dx]['age'] = []
                if loc_list[i] == "N/A":
                    loc_list[i] = "Not Reported"
                if loc_list[i] not in temp_dx[cur_dx]['loc']:
                    temp_dx[cur_dx]['loc'].append(loc_list[i])
                if cdx_list[i] not in temp_dx[cur_dx]['dx']:
                    temp_dx[cur_dx]['dx'].append(cdx_list[i])
                # age = 999
                # if age_list[i] != 'NULL' and age_list[i] != 'None':
                #     age = math.floor(float(age_list[i]) / 365.25)
                # pdb.set_trace()
                if age_list[i] not in temp_dx[cur_dx]['age']:
                    temp_dx[cur_dx]['age'].append(age_list[i])
            else:
                sys.stderr.write('WARN: biospecimen ' + bs_id + ' with dx ' + cdx_list[i] + ' is invalid, skipping!\n')

        for cur_dx in temp_dx:

            if cur_dx not in master_dict:
                master_dict[cur_dx] = {}
                blacklist[cur_dx] = {}
            ages = temp_dx[cur_dx]['age']

            if pt_id not in master_dict[cur_dx]:
                master_dict[cur_dx][pt_id] = {}
                cur_pt = master_dict[cur_dx][pt_id]
                # set to a crazy number, unlikely to see a 100 year old pediatric patient, will adjust this at the end
                cur_pt['age'] = 36525
                cur_pt['gender'] = tumor_object['gender']
                cur_pt['ethnicity'] = tumor_object['ethnicity']
                cur_pt['race'] = tumor_object['race']
                cur_pt['external_id'] = tumor_object['external_id']
                cur_pt['tumor_site'] = ';'.join(temp_dx[cur_dx]['loc'])
                cur_pt['os_age_mos'] = ''
                if tumor_object.get('outcome_age_at_event_days') is not None:
                    try:
                        int(tumor_object['outcome_age_at_event_days'])
                        cur_pt['os_age_mos'] = tumor_object['outcome_age_at_event_days']
                    except Exception as e:
                        print(tumor_object)
                        sys.stderr.write(str(e) + "\nSurvival status age " + tumor_object['outcome_age_at_event_days'] + " not a number.  Setting blank\n")
                        cur_pt['os_age_mos'] = ''

                v_status = ''
                if tumor_object.get('outcome_vital_status') is not None:
                    if tumor_object['outcome_vital_status'] in v_status_dict:
                        v_status = v_status_dict[tumor_object['outcome_vital_status']]
                cur_pt['vital_status'] = v_status
                cur_pt['samples'] = {}
            cur_pt = master_dict[cur_dx][pt_id]
            # calc min age that is not null
            n_flag = 1
            # will convert age to years at end, after earliest dx age determined
            for age in ages:
                try:
                    int(age)
                    n_flag = 0
                    if int(age) < int(cur_pt['age']):
                        cur_pt['age'] = int(age)
                except Exception as e:
                    sys.stderr.write(str(e) + "\nAge could " + age + " not be converted, skipping.\n")
            if samp_id not in cur_pt['samples']:
                cur_pt['samples'][samp_id] = {}
            cur_samp = master_dict[cur_dx][pt_id]['samples'][samp_id]
            if ana_type in cur_samp:
                if samp_id in blacklist[cur_dx] and ana_type != blacklist[cur_dx][samp_id]:
                    blacklist[cur_dx][samp_id] = 'ALL'
                else:
                    blacklist[cur_dx][samp_id] = ana_type
                sys.stderr.write('WARN: Two or more biospecimens of the same analyte type ' + ana_type
                                 + ' share sample ID ' + samp_id + ', seen in analyte types '
                                 + blacklist[cur_dx][samp_id] + ' so far, skipping!\n')
            else:
                cur_samp[ana_type] = {}
                cur_samp[ana_type]['specimen_id'] = bs_id
                cur_samp[ana_type]['tumor_site'] = ';'.join(temp_dx[cur_dx]['loc'])
                cur_samp[ana_type]['cancer_type'] = ';'.join(temp_dx[cur_dx]['dx'])
                tumor_type = tumor_object['source_text_tumor_descriptor']
                cur_samp[ana_type]['tumor_type'] = tumor_type
                cur_samp[ana_type]['sample_type'] = samp_type
                if ana_type == 'DNA':
                    norm_spec = dna_pairs[bs_id]
                    norm_samp = norm_samp_id[norm_spec]
                    cur_samp[ana_type]['matched_norm_samp'] = norm_samp
                    cur_samp[ana_type]['matched_norm_spec'] = norm_spec
                else:
                    cur_samp[ana_type]['rsem'] = bs_id # build_samp_id(rna_samp_def, t_header, info)
    
    for dx in master_dict:
        sys.stderr.write('Outputting results for ' + dx + '\n')
        os.mkdir(out_dir + dx)
        out_samp = open(out_dir + dx + '/data_clinical_sample.txt', 'w')
        out_pt = open(out_dir + dx + '/data_clinical_patient.txt', 'w')
        
        out_samp.write(samp_head)
        out_pt.write(pt_head)
        for pt_id in sorted(master_dict[dx]):
            # track DNA/RNA pairs, adjust sample ID where needed for matched case, unmatched aliquot
            temp = {}
            f = 0
            cur_pt = master_dict[dx][pt_id]
            for samp_id in sorted(cur_pt['samples']):
                f2 = 0
                cur_samp = cur_pt['samples'][samp_id]
                norm_samp = ''
                norm_spec = ''
                cdx = ''
                loc = ''
                tumor_type = ''
                samp_type = ''
                bs_ids = []
                if 'DNA' in cur_samp and (samp_id not in blacklist[dx] or (samp_id in blacklist[dx] and blacklist[dx][samp_id] == 'RNA')):
                    f += 1
                    f2 += 1
                    temp[samp_id] = {}
                    temp[samp_id]['ct'] = f2
                    temp[samp_id]['DNA'] = 1
                    bs_ids.append(cur_samp['DNA']['specimen_id'])
                    norm_samp = cur_samp['DNA']['matched_norm_samp']
                    norm_spec = cur_samp['DNA']['matched_norm_spec']
                    cdx = cur_samp['DNA']['cancer_type']
                    loc = cur_samp['DNA']['tumor_site']
                    tumor_type = cur_samp['DNA']['tumor_type']
                    samp_type = cur_samp['DNA']['sample_type']
                if 'RNA' in cur_samp and (samp_id not in blacklist[dx] or (samp_id in blacklist[dx] and blacklist[dx][samp_id] == 'DNA')):
                    if 'DNA' not in cur_samp or (samp_id in blacklist[dx] and blacklist[dx][samp_id] == 'DNA'):
                        cdx = cur_samp['RNA']['cancer_type']
                        loc = cur_samp['RNA']['tumor_site']
                        tumor_type = cur_samp['RNA']['tumor_type']
                        samp_type = cur_samp['RNA']['sample_type']
                        temp[samp_id] = {}
                    f += 1
                    f2 += 1
                    temp[samp_id]['ct'] = f2
                    temp[samp_id]['RNA'] = 1
                    bs_ids.append(cur_samp['RNA']['specimen_id'])

                if f2 > 0:
                    temp[samp_id]['entry'] = '\t'.join((pt_id, samp_id, ';'.join(bs_ids), cdx, cdx, loc, tumor_type, samp_type, norm_samp, norm_spec, samp_id)) + '\n'
            check = {}
            for samp_id in temp:
                # if cbttc DNA and RNA is paired, or is not cbttc-related or has no RNA print as-is, if not, check id a
                # matched case ID, but unmatched aliqout ID exists
                if temp[samp_id]['ct'] == 2 or config_data['rna_flag'] == 0 \
                        or config_data['rna_samp_id_style'] != 'cbttc_rna_std':
                    out_samp.write(temp[samp_id]['entry'])
                else:
                    check[samp_id] = 0
            # if only one entry in check, no need to look for unmatched pairs, otherwise look
            if len(check) == 1:
                for key in check:
                    out_samp.write(temp[key]['entry'])
            elif len(check) > 1:
                new_ids = {}
                for samp_id in check:
                    parts = samp_id.split('-')
                    new_id = '-'.join(parts[0:3])
                    if re.search('-CL', samp_id):
                        m = re.search('-CL.*', samp_id)
                        new_id += m.group(0)
                    if new_id not in new_ids:
                        new_ids[new_id] = []
                    new_ids[new_id].append(samp_id)
                # now see if new IDs can reveal mismatched aliquots between DNA/RNA
                for new_id in new_ids:
                    if len(new_ids[new_id]) == 1:
                        old_id = new_ids[new_id][0]
                        out_samp.write(temp[old_id]['entry'])
                    else:
                        d = 0
                        r = 0
                        d_id = ''
                        r_id = ''
                        for old_id in new_ids[new_id]:
                            if 'DNA' in temp[old_id]:
                                d += 1
                                d_id = old_id
                            if 'RNA' in temp[old_id]:
                                r += 1
                                r_id = old_id
                        if r == 1 and d == 1:
                            sys.stderr.write('Mismatched aliquots for case ID ' + new_id + ' found, renaming sample\n')
                            d_dict = temp[d_id]['entry'].rstrip('\n').split('\t')
                            r_dict = temp[r_id]['entry'].rstrip('\n').split('\t')
                            new_samp = d_dict
                            new_samp[1] = new_id
                            new_samp[2] = ';'.join((d_dict[2], r_dict[2]))
                            new_samp[10] = d_id + ';' + r_id
                            out_samp.write('\t'.join(new_samp) + '\n')
                        else:
                            for old_id in new_ids[new_id]:
                                out_samp.write(temp[old_id]['entry'])

            if f > 0:
                cur_pt = master_dict[dx][pt_id]
                if cur_pt['age'] == 36525:
                    cur_pt['age'] = ''
                else:
                    age_in_days = cur_pt['age']
                    cur_pt['age'] = str(math.floor(float(age_in_days)/365.25))
                    try: 
                        int(cur_pt['os_age_mos'])
                        diff = int(cur_pt['os_age_mos']) - age_in_days
                        if diff < 0:
                            sys.stderr.write('WARN: OS status occurs before patient dx for ' + pt_id + ' skipping outcome age calc\n')
                            cur_pt['os_age_mos'] = ''
                        elif diff == 0 and cur_pt['vital_status'] == 'LIVING':
                            cur_pt['os_age_mos'] = ''
                        else:
                            cur_pt['os_age_mos'] = str(math.floor(float(diff) / (365.25/12)))
                    except Exception as e:
                        sys.stderr.write(str(e) + "\nSurvival status age " + cur_pt['os_age_mos'] + " not a number.  Setting blank\n")
                        cur_pt['os_age_mos'] = ''
                out_pt.write('\t'.join((pt_id, cur_pt['external_id'], cur_pt['gender'], cur_pt['age'], cur_pt['tumor_site'],
                                cur_pt['race'], cur_pt['ethnicity'], cur_pt['vital_status'], cur_pt['os_age_mos'])) + '\n')
            else:
                sys.stderr.write('WARN: ' + pt_id + ' skipped, all samples were in blacklist\n')
        out_samp.close()
        out_pt.close()

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


bs_attrs = ['external_aliquot_id', 'analyte_type', 'source_text_tissue_type', 'source_text_tumor_descriptor',
            'composition', 'external_sample_id']
pt_attrs = ['external_id', 'gender', 'ethnicity', 'race']
outcome_attrs = ['age_at_event_days', 'vital_status']
dx_attrs = ['source_text_diagnosis', 'source_text_tumor_location', 'age_at_event_days']

x = 1
m = 100
norm_seen = {}
norm_out_fh = []
tum_out_fh = []
v_status_dict = {'Alive': 'LIVING', 'Deceased': 'DECEASED'}
config_data =   {}
out_dir = '/tmp/datasheets/'

# END: STEP 1

if __name__ == '__main__':
    config_file = Path(__file__).parent / "../REFS/data_processing_config.json"
    with open(config_file) as f:
        config_data = json.load(f)

    valid_extensions = []
    rna_ext_dict = config_data['rna_ext_list']
    for key in rna_ext_dict:
        valid_extensions.append(rna_ext_dict[key])
    dna_ext_dict = config_data['dna_ext_list']
    for key in dna_ext_dict:
        valid_extensions.append(dna_ext_dict[key])

    try:
        os.makedirs(out_dir)
    except:
        sys.stderr.write(out_dir + ' already exists.\n')

    pt_head = '\n'.join(config_data['ds_pt_desc'])

    # IMPORTANT! will use external sample id as sample id, and bs id as a specimen id
    samp_head = '\n'.join(config_data['ds_samp_desc'])

    # api = sbg.Api(url='https://cavatica-api.sbgenomics.com/v2', token='296f647e655b4ff2acbc06f92a56b733')

    # project:Project = api.projects.get(id='zhangb1/pnoc003-cbio-data')
    # #TODO: multiple project
    # #projects = api.projects.query(owner='zhangb1', name='PNOC003 Cbio Data')
    # folders: List[File] = list(project.get_files())

    # resources: List[File] =[]

    # for folder in folders:
    #     if folder.is_folder:
    #         # TODO: getting only 5 file fore testing purpose
    #         response = api.files.bulk_get(folder.list_files())
    #         for record in response:
    #             if record.valid:
    #                 name = record.resource.name
    #                 parts = name.split('.')
    #                 ext = ".".join(parts[1:])
    #                 if ext in valid_extensions:
    #                     resources.append(record.resource)
    #                     #print(record.resource.__dict__)
    #             else:
    #                 print(record.error)

    # filteredResources,dna_pairs = get_resources(resources, config_data)

    filteredResources = get_resources_from_cavatica_projects(['zhangb1/pnoc003-cbio-data'], config_data, {})
    dna_pairs = get_tn_pair_data(filteredResources)

    tum_out_fh, norm_out_fh = get_resource_information(filteredResources, config_data)

    # with concurrent.futures.ThreadPoolExecutor(config_data['threads']) as executor:

    #     results = {executor.submit(mt_query_calls, resource): resource for resource in filteredResources}
    #     for result in concurrent.futures.as_completed(results):
    #         if x % m == 0:
    #             sys.stderr.write('Processed ' + str(x) + ' lines\n')
    #             sys.stderr.flush()
    #         (tum_info, norm_info) = result.result()
    #         if bool(tum_info):
    #             tum_out_fh.append(tum_info)
    #         if bool(norm_info):
    #             norm_out_fh.append(norm_info)
    #         x += 1

    dx_dict = get_diagnosis_cbio_disease_mapping(config_data['dx_tbl_fn'])
    norm_samp_id = get_norm_info(norm_out_fh)
    create_master_dict(tum_out_fh, dx_dict, norm_samp_id)

    
