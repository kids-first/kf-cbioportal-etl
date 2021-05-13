#!/usr/bin/env python3

import argparse
import sys
import os
import math
import json
import re
import pdb
from sample_id_builder_helper import build_samp_id


def get_tn_pair_data(c_tbl, bs_ids_blacklist):
    cav_fh = open(c_tbl)
    next(cav_fh)
    pairs = {}
    for line in cav_fh:
        info = line.rstrip('\n').split('\t')
        bs_id = info[0]
        if info[4] == 'RNA':
            rna_data[bs_id] = 0
        else:
            if bs_id not in pairs:
                pairs[bs_id] = info[1]
            else:
                bs_ids_blacklist[bs_id] = 'Double norm'
                sys.stderr.write('WARN: tumor bs id ' + bs_id + ' already associated with a normal sample.  Skipping!\n')
    cav_fh.close()
    return pairs


def get_norm_info(n_tbl):
    dna_samp_def = config_data['dna_norm_id_style']
    norm_fn = n_tbl
    norm_fh = open(norm_fn)
    n_head = next(norm_fh)
    n_header = n_head.rstrip('\n').split('\t')
    n_dict = {}
    bidx = n_header.index('BS_ID')
    eidx = n_header.index('external_sample_id')
    for line in norm_fh:
        info = line.rstrip('\n').split('\t')
        try:
            nsamp = build_samp_id(dna_samp_def, n_header, info)
            n_dict[info[bidx]] = nsamp
        except Exception as e:
            sys.stderr.write(str(e) + ' no normal data skip ' + line)
    norm_fh.close()
    return n_dict


def create_dx_index(d_tbl):
    dx_fh = open(d_tbl)
    next(dx_fh)
    d_dict= {}
    for line in dx_fh:
        (cbttc, cbio_short, cbio_long) = line.rstrip('\n').split('\t')
        d_dict[cbttc] = cbio_short
    dx_fh.close()
    return d_dict


def create_master_dict(t_tbl, dx_dict, norm_samp_id, blacklist, bs_ids_blacklist):
    master_dict = {}
    dna_samp_def = config_data['dna_samp_id_style']
    rna_samp_def = config_data['rna_samp_id_style']
    tum_fh = open(t_tbl)
    # establish indexes of field locations
    t_head = next(tum_fh)
    t_header = t_head.rstrip('\n').split('\t')
    pidx = t_header.index('PT_ID')
    bidx = t_header.index('BS_ID')
    e_pidx = t_header.index('external_id')
    ana_idx = t_header.index('analyte_type')
    sty_idx = t_header.index('composition')
    tty_idx = t_header.index('source_text_tumor_descriptor')
    age_idx = t_header.index('age_at_event_days')
    dx_idx = t_header.index('source_text_diagnosis')
    loc_idx = t_header.index('source_text_tumor_location')
    g_idx = t_header.index('gender')
    eth_idx = t_header.index('ethnicity')
    r_idx = t_header.index('race')
    outcome_age_idx = t_header.index('outcome_age_at_event_days')
    vit_idx = t_header.index('vital_status')
    # collate tum info with DNA pair info and norm info by dx
    sys.stderr.write('Collating information for data sheet\n')
    for line in tum_fh:
        info = line.rstrip('\n').split('\t')
        bs_id = info[bidx]
        if bs_id in bs_ids_blacklist or len(info) < 4:
            continue
        pt_id = info[pidx]
        ana_type = info[ana_idx]
        samp_id = build_samp_id(dna_samp_def, t_header, info)
        samp_type = info[sty_idx]
        if samp_type == 'Derived Cell Line' or samp_type == 'Tissue Cell Culture':
            samp_id += '-CL'
            if args.cl_supp is not None and bs_id in cl_supp:
                samp_id += '-' + cl_supp[bs_id]
        cdx_list = info[dx_idx].split(';')
        loc_list = info[loc_idx].split(';')
        age_list = info[age_idx].split(';')
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
                bs_ids_blacklist[bs_id] = 'No dx\n'
                sys.stderr.write('WARN: biospecimen ' + bs_id + ' with dx\t' + cdx_list[i] + '\tis invalid, skipping!\n')

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
                cur_pt['gender'] = info[g_idx]
                cur_pt['ethnicity'] = info[eth_idx]
                cur_pt['race'] = info[r_idx]
                cur_pt['external_id'] = info[e_pidx]
                cur_pt['tumor_site'] = ';'.join(temp_dx[cur_dx]['loc'])
                cur_pt['os_age_mos'] = ''
                try:
                    int(info[outcome_age_idx])
                    cur_pt['os_age_mos'] = info[outcome_age_idx]
                except Exception as e:
                    sys.stderr.write(str(e) + "\nSurvival status age " + info[outcome_age_idx] + " not a number.  Setting blank\n")
                    cur_pt['os_age_mos'] = ''

                v_status = ''
                if info[vit_idx] in v_status_dict:
                    v_status = v_status_dict[info[vit_idx]]
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
                tumor_type = info[tty_idx]
                cur_samp[ana_type]['tumor_type'] = tumor_type
                cur_samp[ana_type]['sample_type'] = samp_type
                if ana_type == 'DNA':
                    norm_spec = dna_pairs[bs_id]
                    norm_samp = norm_samp_id[norm_spec]
                    cur_samp[ana_type]['matched_norm_samp'] = norm_samp
                    cur_samp[ana_type]['matched_norm_spec'] = norm_spec
                else:
                    cur_samp[ana_type]['rsem'] = bs_id # build_samp_id(rna_samp_def, t_header, info)
    return master_dict, blacklist, bs_ids_blacklist


parser = argparse.ArgumentParser(description='Convert metadata info in data patient and sample sheets for cbio portal')
parser.add_argument('-c', '--cavatica', action='store', dest='cav',
                    help='file with task info from cavatica (see step 1)')
parser.add_argument('-t', '--tumor', action='store', dest='tumor',
                    help='file with tumor metadata (see step 2)')
parser.add_argument('-n', '--normal', action='store', dest='normal', help='file with normal metadata (see step 2)')
parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types and '
                                                                               'data locations')
parser.add_argument('-s', '--cell-line-supplement', action='store', dest='cl_supp', help='supplemental file with cell '
                                                                        'line meta data - bs id<tab>type.  optional')

args = parser.parse_args()
# config_data dict has most customizable options from json config file
with open(args.config_file) as f:
    config_data = json.load(f)
rna_data = {}
out_dir = 'datasheets/'
try:
    os.mkdir(out_dir)
except:
    sys.stderr.write(out_dir + ' already exists.\n')
blacklist = {}
bs_ids_blacklist = {}
v_status_dict = {'Alive': 'LIVING', 'Deceased': 'DECEASED'}
# gathering DNA somatic bs id pairs and RNA bs ids run on cavatica
cav_fn = args.cav
dna_pairs = get_tn_pair_data(args.cav, bs_ids_blacklist)

# gathering norm sample IDS
norm_samp_id = get_norm_info(args.normal)

# create a cbttc dx-to-cbio index
dx_dict = create_dx_index(config_data['dx_tbl_fn'])
# populate cell line supplemental if not none
cl_supp = {}
if args.cl_supp is not None:
    cl_info = open(args.cl_supp)
    for line in cl_info:
        info = line.rstrip('\n').split('\t')
        cl_supp[info[0]] = info[1]

# has all info necessary to flatten by pt and sample, key hierarchy: dx -> pt_id -> pt_attrs, samples -> sample_id -> samp attrs
(master_dict, blacklist, bs_ids_blacklist) = create_master_dict(args.tumor, dx_dict, norm_samp_id, blacklist, bs_ids_blacklist)

sys.stderr.write('Creating dx-specific data sheets\n')
pt_head = '\n'.join(config_data['ds_pt_desc'])

# IMPORTANT! will use external sample id as sample id, and bs id as a specimen id
samp_head = '\n'.join(config_data['ds_samp_desc'])

# output RNA key
pt_RNA = {}
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
        # pt_RNA = {}
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
                pt_RNA[samp_id] = {}
                pt_RNA[samp_id]['ds'] = samp_id
                pt_RNA[samp_id]['bs'] = cur_samp['RNA']['specimen_id']
                pt_RNA[samp_id]['rsem'] = cur_samp['RNA']['rsem']

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
            # pdb.set_trace()
            new_ids = {}
            for samp_id in check:
                # parts = samp_id.split('-')
                # new_id = '-'.join(parts[0:2])
                parts = samp_id.split('_')
                new_id = parts[0]
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
                        # rna_out.write(r_id + '\t' + new_id + '\t' + pt_RNA[r_id]['bs'] + '\n')
                        pt_RNA[r_id]['ds'] = new_id
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
if config_data['rna_flag'] == 1:
    mappingFile = open('mappingFile.txt', 'w')
    for samp_id in pt_RNA:
        mappingFile.write(pt_RNA[samp_id]['ds'] + '\t' + pt_RNA[samp_id]['rsem'] + '\n')
    mappingFile.close()
sys.stderr.write('Outputting complete, check files\n')

