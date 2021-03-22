#!/usr/bin/env python3

import os
import pandas as pd
from pathlib import Path
from scripts.file_utils import write_case_list

case_lists_meta = {
    "cases_3way_complete": {
        "stable_id": "3way_complete",
        "case_list_name": "Tumor samples with mutatation, CNA and mRNA data",
        "case_list_description": "All tumor samples with mutation, CNA, and mRNA data",
        "case_list_category": "all_cases_with_mutation_and_cna_and_mrna_data"
    },
    "cases_all": {
        "stable_id": "all",
        "case_list_name": "All Tumors",
        "case_list_description": "All tumor samples",
        "case_list_category": "all_cases_in_study"
    },
    "cases_cnaseq": {
        "stable_id": "cnaseq",
        "case_list_name": "Tumor samples with mutatation and CNA data",
        "case_list_description": "All tumor samples with mutation and CNA data",
        "case_list_category": "all_cases_with_mutation_and_cna_data"
    },
    "cases_cna": {
        "stable_id": "cna",
        "case_list_name": "Tumor Samples with CNA data",
        "case_list_description": "All tumors with CNA data",
        "case_list_category": "all_cases_with_cna_data"
    },
    "cases_rna_seq_v2_mrna": {
        "stable_id": "rna_seq_v2_mrna",
        "case_list_name": "Tumor Samples with mRNA data (RNA Seq V2)",
        "case_list_description": "All samples with mRNA expression data",
        "case_list_category": "all_cases_with_mrna_rnaseq_data"
    },
    "cases_sequenced": {
        "stable_id": "sequenced",
        "case_list_name": "Tumor samples with mutations",
        "case_list_description": "All tumor samples with mutation data",
        "case_list_category": "all_cases_with_mutation_data"
    },
    "cases_sv": {
        "stable_id": "sv",
        "case_list_name": "Tumor samples with fusions",
        "case_list_description": "All tumor samples with fusion data",
        "case_list_category": "all_cases_with_sv_data"
    }
}


def add_case_list_files(study_id, study_directory):
    case_dir = study_directory + "/case_lists/"
    Path(case_dir).mkdir(exist_ok=True)
    mutations_samples = []
    cna_samples = []

    clinical_sample_fname = study_directory + '/data_clinical_sample.txt'
    if os.path.exists(clinical_sample_fname):
        clinical_sample_df = pd.read_csv(clinical_sample_fname, sep="\t", skiprows=4)
        all_list = clinical_sample_df['SAMPLE_ID'].unique()
        write_case_list(case_dir + 'cases_all.txt', case_lists_meta['cases_all'], study_id, all_list)

    muts_fname = study_directory + '/data_mutations_extended.txt'
    if os.path.exists(muts_fname):
        maf_df = pd.read_csv(muts_fname, sep="\t")
        mutations_samples = maf_df['Tumor_Sample_Barcode'].unique()
        write_case_list(case_dir + 'cases_sequenced.txt', case_lists_meta['cases_sequenced'], study_id,
                        mutations_samples)

    cna_fname = study_directory + '/data_CNA.txt'
    if os.path.exists(cna_fname):
        cna_file = open(cna_fname)
        head = next(cna_file)
        cna_file.close()
        # assumes header is Hugo_symbols\tsample_name1\tsamplename2 etc, if entrez ID, will need to change!
        cna_samples = head.rstrip('\n').split('\t')[1:]
        write_case_list(case_dir + 'cases_cna.txt', case_lists_meta['cases_cna'], study_id, cna_samples)

        if len(mutations_samples) > 0:
            mutation_and_cna_samples = list(set(mutations_samples) & set(cna_samples))
            if len(mutation_and_cna_samples) > 0:
                write_case_list(case_dir + 'cases_cnaseq.txt', case_lists_meta['cases_cnaseq'], study_id,
                                mutation_and_cna_samples)

    rna_fname = study_directory + '/data_rna_seq_v2_mrna_median_Zscores.txt'
    if os.path.exists(rna_fname):
        rna_file = open(rna_fname)
        head = next(rna_file)
        rna_file.close()
        # assumes header is Hugo_symbols\tsample_name1\tsamplename2 etc, if entrez ID, will need to change!
        rna_list = head.rstrip('\n').split('\t')[1:]
        write_case_list(case_dir + 'cases_rna_seq_v2_mrna.txt', case_lists_meta['cases_rna_seq_v2_mrna'], study_id,
                        rna_list)

        if len(mutations_samples) > 0 and len(cna_samples) > 0:
            three_way = list(set(mutations_samples) & set(cna_samples) & set(rna_list))
            if len(three_way) > 0:
                write_case_list(case_dir + 'cases_3way_complete.txt', case_lists_meta['cases_3way_complete'], study_id,
                                three_way)

    fusion_fname = study_directory + '/data_fusions.txt'
    if os.path.exists(fusion_fname):
        fusion_df = pd.read_csv(fusion_fname, sep="\t")
        fusion_list = fusion_df['Tumor_Sample_Barcode'].unique()
        write_case_list(case_dir + 'cases_sv.txt', case_lists_meta['cases_sv'], study_id, fusion_list)
