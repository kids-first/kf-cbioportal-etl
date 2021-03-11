#!/usr/bin/env python3
import sys
import os
import subprocess
import sevenbridges as sbg
import numpy as np
from .process_fusion_data import add_fusion_file
from .process_rsem_data import add_rsem_file
from .process_cnv_data import add_cnv_file
from .process_maf_data import add_maf_file
from .add_case_list_files import add_case_list_files

def add_meta_study_file(meta_data):
    flist = subprocess.check_output('find '+ config_data['datasets'] +' -mindepth 1 -maxdepth 1 -type d', shell=True)

    ds_list = flist.decode().split('\n')
    if ds_list[-1] == '':
        ds_list.pop()

    for dataset_dir in ds_list:
        meta_study = open(dataset_dir + "meta_study.txt", "w")
        meta_study.write("type_of_cancer: " + study + "\n")
        # default is dx value, otherwise make it "cancer_study_identifier" value from config, in not empty
        canc_study_id = study
        if meta_data["cancer_study_identifier"] != "":
            canc_study_id = meta_data["cancer_study_identifier"]
        canc_study_id += meta_data['dir_suffix']
        meta_study.write("cancer_study_identifier: " + canc_study_id + "\n")
        name = dx_dict[study]
        if meta_data["name_append"] != "":
            name += " " + meta_data["name_append"]
        meta_study.write("name: " + name + "\n")
        desc = meta_data["description"][0]
        if len(meta_data["description"]) > 1:
            desc += " " + dx_dict[study] + " " + meta_data["description"][1]
        meta_study.write("description: " + desc + "\nshort_name: " + canc_study_id + "\ngroups: " + meta_data["groups"] + "\n")
        meta_study.close()

def process_ds_new(config_data, resourceSet):
    flist = subprocess.check_output('find '+ config_data['datasets'] +' -name data_clinical_sample.txt', shell=True)

    ds_list = flist.decode().split('\n')
    if ds_list[-1] == '':
        ds_list.pop()

    sbg_api_client = sbg.Api(url='https://cavatica-api.sbgenomics.com/v2', token='296f647e655b4ff2acbc06f92a56b733')
    rsem_files_by_project = {}

    for dpath in ds_list:
        try:
            # last item in find will likely be empty after split
            x = 0
            parts = dpath.split('/')
            cbio_proj = parts[-2]
            # project/disease name should be name of directory hosting datasheet
            sys.stderr.write('Processing ' + cbio_proj + ' project' + '\n')
            cur_ds = open(dpath)
            for i in range(0, 4, 1):
                next(cur_ds)
            ds_head = next(cur_ds)
            ds_header = ds_head.rstrip('\n').split('\t')
            sa_idx = ds_header.index('SAMPLE_ID')
            sp_idx = ds_header.index('SPECIMEN_ID')
            ns_idx = ds_header.index('MATCHED_NORMAL_SPECIMEN_ID')
            studyResources = []
            for entry in cur_ds:
                meta = entry.rstrip('\n').split('\t')
                check = meta[sp_idx].split(';')
                for bs_id in check:
                    # can't tell if RNA or DNA from datasheet, but DNA will be tum + norm bs id, RNA tum + NA
                    key = bs_id + meta[ns_idx]
                    if key in resourceSet:
                        resource = resourceSet[key];
                        resource['SAMPLE_ID'] = meta[sa_idx]
                        studyResources.append(resource)
                    elif bs_id in resourceSet:
                        resource = resourceSet[bs_id];
                        resource['SAMPLE_ID'] = meta[sa_idx]
                        studyResources.append(resource)
                    else:
                        sys.stderr.write('WARNING: ' + key + ' nor ' + bs_id + ' found in data sheet entry\n' + entry + 'Check inputs!\n')

            fileSet = {}
            for studyResource in studyResources:
                for resource in studyResource['resources']:
                    resource.metadata['sample_id'] = studyResource['SAMPLE_ID']
                    name = resource.name
                    name_parts = name.split('.')
                    ext = ".".join(name_parts[1:])

                    if ext in fileSet:
                        fileSet[ext].append(resource)
                    else:
                        fileSet[ext] = [resource]

            maf_extension = config_data['dna_ext_list']['mutation']
            cnv_extension = config_data['dna_ext_list']['copy_number']
            rsem_extension = config_data['rna_ext_list']['expression']
            fusion_extension = config_data['rna_ext_list']['fusion']
            if maf_extension in fileSet:
                add_maf_file(config_data, cbio_proj, "/".join(parts[:-1]),  fileSet[maf_extension])
            if cnv_extension in fileSet:
                add_cnv_file(config_data, sbg_api_client, cbio_proj, "/".join(parts[:-1]),  fileSet[cnv_extension])
            if rsem_extension in fileSet:
                rsem_files_by_project[cbio_proj] = fileSet[rsem_extension]
            if fusion_extension in fileSet:
                add_fusion_file(config_data, "/".join(parts[:-1]),  fileSet[fusion_extension])

        except Exception as e:
            print(e)
            sys.exit()

    if bool(rsem_files_by_project):
        add_rsem_file(config_data, sbg_api_client, rsem_files_by_project)

    add_case_list_files(config_data)

    add_meta_study_file(config_data)
