#!/usr/bin/env python3
"""
This is a genomics file + clinical data file to cBio package conversion workflow script.
It is designed to work with standard KF somatic workflow outputs as well as DGD outputs.
Clinical files should have been produced ahead of time, while supporting sample ID file manifest, case_meta config json file, and data json config.
If there are fusion files, they also need s table outlining sequencing center location.
See README for prequisite details.
"""
from re import sub
import sys
import os
import argparse
import json
import subprocess
import pdb


parser = argparse.ArgumentParser(description='Download files (if needed), collate genomic files, organize load package.')
parser.add_argument('-t', '--table', action='store', dest='table',
                    help='Table with cbio project, kf bs ids, cbio IDs, and file names')
parser.add_argument('-m', '--manifest', action='store', dest='manifest', help='Download file manifest, if needed')
parser.add_argument('-c', '--cbio-config', action='store', dest='cbio_config', help='cbio case and meta config file')
parser.add_argument('-d', '--data-config', action='store', dest='data_config', help='json config file with data types and '
                                                                               'data locations')
parser.add_argument('-f', '--dgd-status', action='store', dest='dgd_status', help='Flag to determine load will have pbta/kf + dgd(both), kf/pbta only(kf), dgd-only(dgd)',
default='both', const='both', nargs='?', choices=['both', 'kf', 'dgd'])


args = parser.parse_args()
with open(args.data_config) as f:
    config_data = json.load(f)


def process_maf():
    """
    Collate and process pbta/kf style mafs. Call append if also adding dgd data
    """
    maf_dir = config_data['file_loc_defs']['mafs']['kf']
    if args.dgd_status == 'dgd':
        config_data['file_loc_defs']['mafs']['dgd']
    maf_header = config_data['file_loc_defs']['mafs']['header']
    maf_cmd = script_dir + 'maf_merge.py -t ' + args.table + ' -i ' + maf_header + ' -m ' + maf_dir + ' -j ' + args.data_config + ' -f ' + args.dgd_status + ' 2> collate_mafs.log'
    status = 0
    status = subprocess.call(maf_cmd, shell=True)
    if args.dgd_status == 'both':
        status += process_append_dgd_maf()
    return status


def process_append_dgd_maf():
    """
    Append DGD mafs to existing collated kf maf
    """
    maf_header = config_data['file_loc_defs']['mafs']['header']
    in_maf_dir = config_data['file_loc_defs']['mafs']['dgd']
    append_maf = 'merged_mafs/' + cbio_study_id + '.maf'
    append_cmd = script_dir + 'add_dgd_maf_to_pbta.py -i ' + maf_header + ' -m ' + in_maf_dir + ' -t ' + args.table + ' >> ' + append_maf + ' 2> dgd_append_maf.log'
    status = 0
    status = subprocess.call(append_cmd, shell=True)
    return status


def process_cnv(cnv_loc_dict):
    """
    Add gene info to CNV calls, merge into table, and create GISTIC-style output
    """
    cnv_dir = cnv_loc_dict['pval']
    gene_annot_cmd = script_dir + 'cnv_1_genome2gene.py -d ' + cnv_dir + ' -j ' +  args.data_config + ' 2> cnv_gene_annot.log'
    status = subprocess.call(gene_annot_cmd, shell=True)
    info_dir = cnv_loc_dict['info']
    merge_gene_cmd = script_dir + 'cnv_2_merge.py -t ' + args.table + ' -n converted_cnvs -j ' + args.data_config
    if info_dir != "":
        merge_gene_cmd += ' -i ' + info_dir
    merge_gene_cmd += ' 2> merge_cnv_gene.log'
    status = 0
    status += subprocess.call(merge_gene_cmd, shell=True)
    gistic_style_cmd = script_dir + 'cnv_3_gistic_style.py -d merged_cnvs -j ' + args.data_config + ' -t ' + args.table
    if info_dir != "":
        gistic_style_cmd += ' -i ' + info_dir
    gistic_style_cmd += ' 2> merge_cnv_gene.log'
    status += subprocess.call(gistic_style_cmd, shell=True)

    if 'seg' in cnv_loc_dict:
        status += process_seg(cnv_loc_dict)
    
    return status


def process_seg(cnv_loc_dict):
    """
    Collate and process CNV seg files
    """
    seg_dir = cnv_loc_dict['seg']
    merge_seg_cmd = script_dir + 'cnv_merge_seg.py -t ' + args.table + ' -m ' + seg_dir + ' -j ' + args.data_config + ' 2> cnv_merge_seg.log'
    status = subprocess.call(merge_seg_cmd, shell=True)
    return status


def process_rsem():
    """
    Merge rsem results by FPKM, calculate z-scores
    """
    rsem_dir = config_data['file_loc_defs']['rsem']
    merge_rsem_cmd = script_dir + 'rna_merge_rename_expression.py -t ' + args.table + ' -r ' + rsem_dir + ' 2> rna_merge_rename_expression.log'
    status = subprocess.call(merge_rsem_cmd, shell=True)
    return status


def process_kf_fusion():
    """
    Collate and process annoFuse output
    """
    fusion_dir = config_data['file_loc_defs']['fusion']
    sq_file = config_data['file_loc_defs']['fusion_sq_file']
    fusion_cmd = script_dir + 'rna_convert_fusion.py -t ' + args.table + ' -f ' + fusion_dir + ' -m annofuse -s ' + sq_file + ' 2> rna_convert_fusion.log'
    status = subprocess.call(fusion_cmd, shell=True)
    return status


def process_dgd_fusion():
    """
    Collate process DGD fusion output
    """
    fusion_dir = config_data['file_loc_defs']['dgd_fusion']
    dgd_fusion_cmd = script_dir + 'add_dgd_fusion.py -t ' + args.table + ' -f ' + fusion_dir
    if args.dgd_status == 'both':
        append_fusion = 'merged_fusion/' + cbio_study_id + '.fusions.txt'
        dgd_fusion_cmd += ' -a >> ' + append_fusion
    dgd_fusion_cmd += ' 2> add_dgd_fusion.log'
    status = subprocess.call(dgd_fusion_cmd, shell=True)
    return status

# This part is commented out for now as there is no good solution yet for files stored in different account - need to be able to run chopaws to get key
# if(args.manifest != None):
#     sys.stderr.write('Download manifest given. Downloading files\n')
#     dl_file_cmd = config_data['script_dir'] + 'get_files_from_manifest.py -m ' + args.manifest + ' -f ' + join(config_data['dl_file_type_list'] + ' -p saml')
#     subprocess.call(dl_file_cmd, shell=True)

with open(args.cbio_config) as f:
    config_meta_case = json.load(f)

cbio_study_id = config_meta_case["study"]["cancer_study_identifier"]
# iterate through config file - file should only have keys related to data to be loaded
script_dir = config_data['script_dir']
for key in config_meta_case:
    if key.startswith('merged_'):
        data_type = '_'.join(key.split('_')[1:])
        if data_type == 'mafs':
            exit_status = process_maf()
            if exit_status:
                sys.stderr.write('Something whent wrong while processing the mafs! \
                and/or dgd_append_maf Check collate_mafs.log for more info\n')
                exit(1)
        elif data_type == 'cnvs':
            exit_status =process_cnv(config_data['file_loc_defs']['cnvs'])
            if exit_status:
                sys.stderr.write('Something whent wrong while processing the cnv data! Check cnv_* logs for more info\n')
                exit(1)
        elif data_type == 'rsem':
            exit_status = process_rsem()
            if exit_status:
                sys.stderr.write('Something whent wrong while processing the rna expression! Check rna_merge_rename_expression.log for more info\n')
                exit(1)
        elif data_type == 'fusion':
            # Status both works for...both, only when one is specificaly picked should one not be run
            if args.dgd_status != 'dgd':
                exit_status = process_kf_fusion()
                if exit_status:
                    sys.stderr.write('Something whent wrong while processing the kf fusions! Check rna_convert_fusion.log for more info\n')
                    exit(1)
            if args.dgd_status != 'kf':
                exit_status = process_dgd_fusion()
                if exit_status:
                    sys.stderr.write('Something whent wrong while processing the dgd fusions! Check add_dgd_fusion.log for more info\n')
                    exit(1)

# Run final package builder script
pck_cmd = config_data['script_dir'] + 'organize_upload_packages.py -o processed -c ' + args.cbio_config + ' 2> load_package_create.log'
exit_status = subprocess.call(pck_cmd, shell=True)
if exit_status:
    sys.stderr.write('Something went wrong while creating the load package. Check load_package_create.log for more info\n')

# Run cbioportal data validator
validate = config_data['cbioportal_validator'] + ' -s processed/' + cbio_study_id + ' -n -v 2> validator.errs > validator.out'
exit_status = subprocess.call(validate, shell=True)
if exit_status:
    sys.stderr.write('Validator quit with status ' + str(exit_status) + '. Check validator.errs for more info\n')
