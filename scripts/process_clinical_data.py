#!/usr/bin/env python3

import sys
from pathlib import Path
import math
import re

from scripts.sample_id_builder_helper import format_smaple_id
from scripts.kf_dataservice_resources import get_tumor_resources, query_dataservice_bs_id

v_status_dict = {'Alive': 'LIVING', 'Deceased': 'DECEASED'}


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


def add_clinical_data_files(config_data, study_id, resources, dx_dict, cl_supp=None):
    sys.stderr.write('Processing for clinical data' + study_id + ' project' + '\n')
    tumor_objects = get_tumor_resources(resources, config_data)
    tumor_bs_mapped_normal_sample = get_tumor_bs_mapped_normal_sample(config_data, resources)

    pt_head = '\n'.join(config_data['ds_pt_desc'])
    # IMPORTANT! will use external sample id as sample id, and bs id as a specimen id
    samp_head = '\n'.join(config_data['ds_samp_desc'])

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
        if samp_type == 'Derived Cell Line':
            samp_id += '-CL'
            if cl_supp is not None and bs_id in cl_supp:
                samp_id += '-' + cl_supp[bs_id]

        if pt_id not in master_dict:
            master_dict[pt_id] = {}
            cur_pt = master_dict[pt_id]
            cur_pt['gender'] = tumor_object['gender']
            cur_pt['ethnicity'] = tumor_object['ethnicity']
            cur_pt['race'] = tumor_object['race']
            cur_pt['external_id'] = tumor_object['external_id']
            cur_pt['tumor_site'] = ";".join(x['source_text_tumor_location'] for x in tumor_object['diagnosis'])
            cur_pt['os_age_mos'] = ''
            if tumor_object.get('outcome_age_at_event_days') is not None:
                try:
                    cur_pt['os_age_mos'] = int(tumor_object['outcome_age_at_event_days'])
                except Exception as e:
                    sys.stderr.write(str(e) + "\nSurvival status age " + tumor_object[
                        'outcome_age_at_event_days'] + " not a number.  Setting blank\n")
                    cur_pt['os_age_mos'] = ''

            v_status = ''
            if tumor_object.get('outcome_vital_status') is not None:
                if tumor_object['outcome_vital_status'] in v_status_dict:
                    v_status = v_status_dict[tumor_object['outcome_vital_status']]
            cur_pt['vital_status'] = v_status
            cur_pt['samples'] = {}

        cur_pt = master_dict[pt_id]
        # calc min age that is not null
        # will convert age to years at end, after earliest dx age determined
        cur_pt['age'] = min(
            x['age_at_event_days'] for x in tumor_object['diagnosis'] if isinstance(x['age_at_event_days'], int))
        if samp_id not in cur_pt['samples']:
            cur_pt['samples'][samp_id] = {}
        cur_samp = master_dict[pt_id]['samples'][samp_id]
        if ana_type in cur_samp:
            if samp_id in blacklist and ana_type != blacklist[samp_id]:
                blacklist[samp_id] = 'ALL'
            else:
                blacklist[samp_id] = ana_type
            sys.stderr.write('WARN: Two or more biospecimens of the same analyte type ' + ana_type
                             + ' share sample ID ' + samp_id + ', seen in analyte types '
                             + blacklist[samp_id] + ' so far, skipping!\n')
        else:
            cur_samp[ana_type] = {}
            cur_samp[ana_type]['specimen_id'] = bs_id
            cur_samp[ana_type]['tumor_site'] = ";".join(
                x['source_text_tumor_location'] for x in tumor_object['diagnosis'])
            cur_samp[ana_type]['cancer_type'] = ";".join(x['source_text_diagnosis'] for x in tumor_object['diagnosis'])
            cur_samp[ana_type]['tumor_type'] = tumor_object['source_text_tumor_descriptor']
            cur_samp[ana_type]['sample_type'] = samp_type
            if ana_type == 'DNA':
                cur_samp[ana_type]['matched_norm_samp'] = tumor_bs_mapped_normal_sample[bs_id]['sample_id']
                cur_samp[ana_type]['matched_norm_spec'] = tumor_bs_mapped_normal_sample[bs_id]['specimen_id']
            else:
                cur_samp[ana_type]['rsem'] = bs_id  # build_samp_id(rna_samp_def, t_header, info)

    Path(config_data['datasets'] + study_id).mkdir(exist_ok=True)
    out_samp = open(config_data['datasets'] + study_id + '/data_clinical_sample.txt', 'w')
    out_pt = open(config_data['datasets'] + study_id + '/data_clinical_patient.txt', 'w')

    out_samp.write(samp_head)
    out_pt.write(pt_head)
    for pt_id in sorted(master_dict):
        # track DNA/RNA pairs, adjust sample ID where needed for matched case, unmatched aliquot
        temp = {}
        f = 0
        cur_pt = master_dict[pt_id]
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
            if 'DNA' in cur_samp and (samp_id not in blacklist or blacklist[samp_id] == 'RNA'):
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
            if 'RNA' in cur_samp and (samp_id not in blacklist or blacklist[samp_id] == 'DNA'):
                if 'DNA' not in cur_samp or (samp_id in blacklist and blacklist[samp_id] == 'DNA'):
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
                temp[samp_id]['entry'] = '\t'.join((pt_id, samp_id, ';'.join(bs_ids), cdx, cdx, loc, tumor_type,
                                                    samp_type, norm_samp, norm_spec, samp_id)) + '\n'
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
            cur_pt = master_dict[pt_id]
            if cur_pt['age'] == 36525:
                cur_pt['age'] = ''
            else:
                age_in_days = cur_pt['age']
                cur_pt['age'] = str(math.floor(float(age_in_days) / 365.25))
                try:
                    int(cur_pt['os_age_mos'])
                    diff = int(cur_pt['os_age_mos']) - age_in_days
                    if diff < 0:
                        sys.stderr.write(
                            'WARN: OS status occurs before patient dx for ' + pt_id + ' skipping outcome age calc\n')
                        cur_pt['os_age_mos'] = ''
                    elif diff == 0 and cur_pt['vital_status'] == 'LIVING':
                        cur_pt['os_age_mos'] = ''
                    else:
                        cur_pt['os_age_mos'] = str(math.floor(float(diff) / (365.25 / 12)))
                except Exception as e:
                    sys.stderr.write(str(e) + "\nSurvival status age " + cur_pt[
                        'os_age_mos'] + " not a number.  Setting blank\n")
                    cur_pt['os_age_mos'] = ''
            out_pt.write(
                '\t'.join((pt_id, cur_pt['external_id'], cur_pt['gender'], cur_pt['age'], cur_pt['tumor_site'],
                           cur_pt['race'], cur_pt['ethnicity'], cur_pt['vital_status'],
                           cur_pt['os_age_mos'])) + '\n')
        else:
            sys.stderr.write('WARN: ' + pt_id + ' skipped, all samples were in blacklist\n')
    out_samp.close()
    out_pt.close()
