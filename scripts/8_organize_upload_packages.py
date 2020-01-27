#!/usr/bin/env python3

import sys
import argparse
import os
import concurrent.futures
import re
import json
import subprocess

parser = argparse.ArgumentParser(description='Create cases lists, meta files, and organize data for cbio upload.'
                                ' It is assumed you are at the dir level of all input data files')
parser.add_argument('-o', '--output_dir', action='store', dest='out_dir', help='output directory name')
parser.add_argument('-d', '--dx-file', action='store', dest='dx_file', help='Dx index file with cbio short and long names')
parser.add_argument('-c', '--config', action='store', dest='config_file', help='json config file with meta information; see REFS/case_meta_config.json example')
args = parser.parse_args()

def process_meta_study(meta_data, output_dir, study):
    meta_study = open(output_dir + "meta_study.txt", "w")
    meta_study.write("type_of_cancer: " + study + "\n")
    meta_study.write("cancer_study_identifier: " + study + "\n")
    # default is dx value, otherwise make it "cancer_study_identifier" value from config, in not empty
    canc_study_id = study
    if meta_data["cancer_study_identifier"] != "":
        canc_study_id = meta_data["cancer_study_identifier"]
    canc_study_id += + meta_data['dir_suffix']
    meta_study.write("canc_study_id: " + canc_study_id + "\n")
    name = dx_dict[study]
    if meta_data["name_append"] != "":
        name += " " + meta_data["name_append"]
    meta_study.write("name: " + name + "\n")
    desc = meta_data["description"][0]
    if len(meta_data["description"]) > 1:
        desc += " " + name + " " + meta_data["description"][1]
    meta_study.write("description: " + desc + "\nshort_name: " + canc_study_id + "\ngroups: " + meta_study["groups\n"])
    meta_study.close()
    return canc_study_id


def process_meta_data(meta_data, output_dir, canc_study_id, study):
    for dtype in meta_data["dtypes"]:
        cbio_name = meta_data["dtypes"]["cbio_name"]
        parts = cbio_name.split("_")
        meta_name = "meta_" + "_".join(parts[1])
        meta_data_file = open(output_dir + meta_name, "w")
        meta_data_file.write("cancer_study_identifier: " + canc_study_id + "\n")
        attr_dict = meta_data["dtypes"]["meta_file_attr"]
        for mkey in attr_dict:
            meta_data_file.write(mkey + ": " + attr_dict[mkey] + "\n")
        meta_data_file.close()
        # create data_ links to data
        cmd = "ln -s " + cwd + meta_data["dir"] + "/" + study + "." + meta_data["dtypes"]["ext"] + " " + output_dir + study + "/" + cbio_name
        subprocess.call(cmd, shell=True)


def process_clinical_data(meta_data, output_dir, canc_study_id, study):
    for dtype in meta_data["dtypes"]:
        cbio_name = meta_data["dtypes"]["cbio_name"]
        parts = cbio_name.split("_")
        meta_name = "meta_" + "_".join(parts[1])
        meta_data_file = open(output_dir + meta_name, "w")
        meta_data_file.write("cancer_study_identifier: " + canc_study_id + "\n")
        attr_dict = meta_data["dtypes"]["meta_file_attr"]
        for mkey in attr_dict:
            meta_data_file.write(mkey + ": " + attr_dict[mkey] + "\n")
        meta_data_file.close()
        # create data_ links to data
        cmd = "ln -s " + cwd + meta_data["dir"] + "/" + study + "/" + cbio_name + " " + output_dir + study + "/" + cbio_name
        subprocess.call(cmd, shell=True)


def create_case_lists(data_dict, output_dir, canc_study_id, study):
    try:
        case_dir = out_dir + study + "/case_lists/"
        os.mkdir(case_dir)
    except:
        sys.stderr.write(case_dir + ' already exists.\n')

    muts_list = []
    cna_list = []
    rna_list = []
    fusion_list = []

    muts_fname = output_dir + study + "/" + config_data["merged_mafs"]["dtypes"]["mutation"]["cbio_name"]
    muts_file = open(muts_fname)
    head = next(muts_file)
    head = next(muts_file)
    header = head.rstrip('\n').split('\t')
    s_idx = header.index('Tumor_Sample_Barcode')
    for line in muts_file:
        data = line.rstrip('\n').spliot('\t')
        muts_list.append(data[s_idx])
    muts_file.close()
    muts_list = [*{*muts_list}]
    if data_dict["merged_cnvs"] == 1:
        cna_fname = output_dir + study + "/" + config_data["merged_cnvs"]["dtypes"]["linear"]["cbio_name"]
        cna_file = open(cna_fname)
        head = next(cna_file)
        # assumes header is Hugo_symbols\tsample_name1\tsamplename2 etc, if entrez ID, will need to change!
        cna_list = head.rstrip('\n').split('\t')[1:]
        cna_file.close()
    if data_dict["merged_rsem"] == 1:
        # assumes header is Hugo_symbols\tsample_name1\tsamplename2 etc, if entrez ID, will need to change!
        rna_fname = output_dir + study + "/" + config_data["merged_rsem"]["dtypes"]["counts"]["cbio_name"]
        rna_file = open(rna_fname)
        head = next(rna_file)
        rna_list = head.rstrip('\n').split('\t')[1:]
    if data_dict["merged_fusion"] == 1:
        fusion_fname = output_dir + study + "/" + config_data["merged_fusion"]["dtypes"]["fusion"]["cbio_name"]
        fusion_file = open(fusion_fname)
        head = next(fusion_file)
        header = head.rstrip('\n').split('\t')
        s_idx = header.index('Tumor_Sample_Barcode')
        for line in fusion_file:
            data = line.rstrip('\n').spliot('\t')
            fusion_list.append(data[s_idx])
        fusion_file.close()
        fusion_list = [*{*fusion_list}]
    




cwd = os.getcwd() + "/"
with open(args.config_file) as f:
    config_data = json.load(f)

dx_file = open(args.dx_file)
dx_dict = {}
head = next(dx_file)
for line in dx_file:
    (ds_dx, cbio_short, cbio_long) = line.rstrip('\n').split('\t')
    dx_dict[cbio_short] = cbio_long
dx_file.close()

out_dir = args.out_dir
if out_dir[-1] != "/":
    out_dir += "/"
try:
    os.mkdir(out_dir)
except:
    sys.stderr.write(out_dir + ' already exists.\n')

for dx in dx_dict:
    cur_dir = out_dir + dx
    if config_data['study']['dir_suffix'] != "":
        cur_dir += config_data['study']['dir_suffix']
    cur_dir += "/"
    try:
        os.mkdir(cur_dir)
    except:
        sys.stderr.write(cur_dir + ' already exists.\n')
    sys.stderr.write("Creating meta study file for " + dx + "\n")
    canc_study_id = process_meta_study(config_data['study'], cur_dir, dx)
    data_keys = {"merged_mafs": 0, "merged_cnvs": 0, "merged_rsem": 0, "merged_fusion": 0}
    for key in data_keys:
        if config_data[key]["dir"] != "":
            data_keys[key] = 1
            process_meta_data(config_data[key], cur_dir, canc_study_id, dx)
            sys.stderr.write("Creating meta data files and links for " + key + "\n")
        else:
            sys.stderr.write("Skipping meta files for " + key + "\n")
    sys.stderr.write("Creating clinical meta sheets and links\n")
    process_clinical_data(config_data["data_sheets"], cur_dir, canc_study_id, dx)