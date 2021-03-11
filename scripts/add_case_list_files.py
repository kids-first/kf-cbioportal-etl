#!/usr/bin/env python3
import sys
import os
import subprocess
import pandas as pd
from pathlib import Path

from .file_utils import write_case_list


def add_case_list_files(config_data):
    flist = subprocess.check_output('find '+ config_data['datasets'] +' -mindepth 1 -maxdepth 1 -type d', shell=True)

    ds_list = flist.decode().split('\n')
    if ds_list[-1] == '':
        ds_list.pop()

    for dpath in ds_list:
        create_case_lists(config_data, dpath)

def create_case_lists(config_data, dataset_dir):
    cbio_proj = dataset_dir.split('/')[-1]
    case_dir = dataset_dir + "/case_lists/"
    Path(case_dir).mkdir(exist_ok=True)

    muts_list = []
    cna_list = []

    clinical_sample_fname = dataset_dir + '/data_clinical_sample.txt'
    if os.path.exists(clinical_sample_fname):
        clinical_sample_df = pd.read_csv(clinical_sample_fname, sep="\t", skiprows=4)
        all_list = clinical_sample_df['SAMPLE_ID'].unique()
        write_case_list(case_dir+'cases_all.txt', config_data['case_lists']['cases_all'], cbio_proj, all_list)

    muts_fname = dataset_dir + '/data_mutations_extended.txt'
    if os.path.exists(muts_fname):
        maf_df = pd.read_csv(muts_fname, sep="\t")
        muts_list = maf_df['Tumor_Sample_Barcode'].unique()
        write_case_list(case_dir+'cases_sequenced.txt', config_data['case_lists']['cases_sequenced'], cbio_proj, muts_list)

    cna_fname = dataset_dir + '/data_CNA.txt'
    if os.path.exists(cna_fname):
        cna_file = open(cna_fname)
        head = next(cna_file)
        cna_file.close()
        # assumes header is Hugo_symbols\tsample_name1\tsamplename2 etc, if entrez ID, will need to change!
        cna_list = head.rstrip('\n').split('\t')[1:]
        write_case_list(case_dir+'cases_cna.txt', config_data['case_lists']['cases_cna'], cbio_proj, cna_list)

        if len(muts_list) > 0:
            muts_plus_cna = list(set(muts_list) & set(cna_list))
            if len(muts_plus_cna) > 0:
                write_case_list(case_dir+'cases_cnaseq.txt', config_data['case_lists']['cases_cnaseq'], cbio_proj, muts_plus_cna)

    rna_fname = dataset_dir + '/data_rna_seq_v2_mrna_median_Zscores.txt'
    if os.path.exists(rna_fname):
        rna_file = open(rna_fname)
        head = next(rna_file)
        rna_file.close()
        # assumes header is Hugo_symbols\tsample_name1\tsamplename2 etc, if entrez ID, will need to change!
        rna_list = head.rstrip('\n').split('\t')[1:]
        write_case_list(case_dir+'cases_rna_seq_v2_mrna.txt', config_data['case_lists']['cases_rna_seq_v2_mrna'], cbio_proj, rna_list)

        if len(muts_list) > 0 and len(cna_list) > 0:
            three_way = list(set(muts_list) & set(cna_list) & set(rna_list))
            if len(three_way) > 0:
                write_case_list(case_dir+'cases_3way_complete.txt', config_data['case_lists']['cases_3way_complete'], cbio_proj, three_way)


    fusion_fname = dataset_dir + '/data_fusions.txt'
    if os.path.exists(fusion_fname):
        fusion_df = pd.read_csv(fusion_fname, sep="\t")
        fusion_list = fusion_df['Tumor_Sample_Barcode'].unique()
        write_case_list(case_dir+'cases_sv.txt', config_data['case_lists']['cases_sv'], cbio_proj, fusion_list)
