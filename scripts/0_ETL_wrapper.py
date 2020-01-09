#!/usr/bin/env python3

import sys
import subprocess
import json
import argparse
import os
from argparse import RawTextHelpFormatter


def check_vars(var_list, step_sub_list):
    flag = 0
    for var in var_list:
        if var not in config_data:
            flag = 1
            sys.stderr.write('ERROR: ' + var + ' not defined in ' + config_file + ' for steps ' + step_sub_list + '\n')
    if flag:
        sys.exit(1)


def config_precheck():
    if 'all' in steps:
        var_list = ['threads', 'cpus']
        check_vars(var_list, 'all')
    if '1' in steps or '4' in steps or '5a' in steps or '5b' in steps or '6' in steps:
        var_list = ['cpus']
        check_vars(var_list, '1, 4, 5a, 5b, 6')
    if '2' in steps:
        var_list = ['threads']
        check_vars(var_list, '2')
    if 'all' in steps or '1' in steps or '3' in steps or '7' in steps:
        var_list = ['rna_flag']
        check_vars(var_list, 'all, 1, 7')
    if 'all' in steps or '1' in steps:
        var_list = ['cav_dna_projects', 'tn_order']
        check_vars(var_list, 'all, 1')
        if config_data['rna_flag']:
            check_vars(['cav_rna_projects'], 'all, 1')
    if 'all' in steps or '1' in steps or '7' in steps:
        var_list = ['cna_flag']
        check_vars(var_list, 'all, 1, 7')
    if 'all' in steps or '3' in steps:
        var_list = ['dna_samp_id_style', 'dx_tbl_fn', 'ds_pt_desc', 'ds_samp_desc']
        if config_data['rna_flag']:
            var_list.append('rna_samp_id_style')
        check_vars(var_list, 'all, 3')
        input_config_files = [config_data['dx_tbl_fn']]
        check_files(input_config_files, 'CONFIG', '3')
    if 'all' in steps or '5a' in steps:
        var_list = ['bedtools', 'cp_only_script', 'bed_genes', 'hugo_tsv']
        check_vars(var_list, 'all, 5a')
        input_config_files = [config_data['cp_only_script'], config_data['bed_genes'],
                              config_data['hugo_tsv']]
        check_files(input_config_files, 'CONFIG', '5a')
    if 'all' in steps or '6' in steps:
        var_list = ['ens_gene_list', 'mapFileLoc']
        check_vars(var_list, 'all, 6')
        input_config_files = [config_data['ens_gene_list'], config_data['mapFileLoc']]
        check_files(input_config_files, 'CONFIG', '6')
    if 'all' in steps or '7' in steps:
        var_list = ['dx_tbl_fn', 'data_dir', 'cancStudyID', 'group', 'study_short_desc', 'study_desc_p1',
                    'study_desc_p2']
        if config_data['rna_flag']:
            var_list.extend(['fpkmFileLoc', 'mapFileLoc', 'genesFileLoc'])
        check_vars(var_list, 'all, 7')
        if config_data['rna_flag']:
            input_config_files = [config_data['mapFileLoc'], config_data['genesFileLoc']]
            check_files(input_config_files, 'CONFIG', '7')
        input_config_dirs = [config_data['data_dir']]
        check_dirs(input_config_dirs, 'CONFIG', '7')


def check_files(file_list, ftype, step):
    flag = 0
    for fn in file_list:
        if fn is None:
            if step[0] != '7' and ftype == 'INPUT':
                sys.stderr.write('ERROR: Missing an input file completely. Try running ' + script_dir + step + ' -h\n')
            else:
                sys.stderr.write('ERROR: No file defined of type ' + ftype + ' for step ' + step + '\n')
            flag = 1
        elif not os.path.isfile(fn):
            sys.stderr.write('ERROR: ' + ftype + ' file ' + fn + ' not found. Check file names and config file\n')
            flag = 1
    if flag == 1:
        sys.stderr.write('Missing files for step(s) ' + step + '.  Exiting.')
        sys.exit(1)


def check_dirs(dir_list, ftype, step):
    flag = 0
    for dn in dir_list:
        if dn is None:
            if step[0] != '7' and ftype == 'INPUT':
                sys.stderr.write('ERROR: Missing an input dir  completely. Try running ' + script_dir + step + ' -h\n')
            else:
                sys.stderr.write('ERROR: No dir defined of type ' + ftype + ' for step ' + step + '\n')
            flag = 1
        elif not os.path.isdir(dn):
            sys.stderr.write('ERROR: ' + ftype + ' directory ' + dn
                             + ' not found. Check directory names and config file\n')
            flag = 1
    if flag == 1:
        sys.stderr.write('Missing directories for step(s) ' + step + '. Exiting.\n')
        sys.exit(1)


parser = argparse.ArgumentParser(description='Wrapper script to run all or some of ETL component scripts.',
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument('-x', '--steps', action='store', dest='steps', help='csv list of steps to execute, valid choices '
        'are: "all", "1", "2", "3", "4", "5a", "5b", "6", "7".  Refers to scripts:\n'
        '1_query_cavatica.py\n'
        '2_query_ds_by_bs_id.py\n'
        '3_create_data_sheets.py\n'
        '4_merge_maf.py\n'
        '5a_cnv_genome2gene.py\n'
        '5b_merge_cnv.py\n'
        '6_merge_rsem.py\n'
        '7_createLoadingFiles.R')
parser.add_argument('-c', '--cavatica', action='store', dest='cav', help='cavatica output file name, relevant steps: '
                                                                         '1, 2, 3, 4, 5b, 6')
parser.add_argument('-p', '--profile', action='store', dest='profile', help='cavatica profile name. requires '
                                        '.sevenbridges/credentials file be present. Relevant steps: 1, 6')
parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types and '
                                                                               'data locations')
parser.add_argument('-u', '--kf-url', action='store', dest='url',
                    help='Kids First data service url, i.e. https://kf-api-dataservice.kidsfirstdrc.org/. '
                         'Relevant step: 2')
parser.add_argument('-t', '--tumor', action='store', dest='tumor',
                    help='file with tumor metadata (see step 2).  Relevant step: 3, not needed if step 2 run')
parser.add_argument('-n', '--normal', action='store', dest='normal', help='file with normal metadata (see step 2)'
                                                                          ' Relevant step: 3, not needed if step 2 run')
parser.add_argument('-s', '--cell-line-supplement', action='store', dest='cl_supp', help='supplemental file with cell '
                                                        'line meta data - bs id<tab>type. Relevant step: 2, OPTIONAL')
parser.add_argument('-i', '--header', action='store', dest='header', help='File with maf header only. Relevant step: 4')
parser.add_argument('-m', '--maf-dir', action='store', dest='maf_dir', help='maf file directory. Relevant step: 4')
parser.add_argument('-d', '--cnv-dir', action='store', dest='cnv_dir', help='cnv as genome coords file directory.'
                                                                            ' Relevant step: 5a')
parser.add_argument('-e', '--cnv-dir-gene', action='store', dest='cnv_dir_gene', help='cnv as gene file directory.'
                                                                          'Relevant step: 5b, not needed if step5a run')
parser.add_argument('-r', '--rsem-dir', action='store', dest='rsem_dir', help='rsem file directory. Relevant step: 6')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

steps = args.steps.split(',')
valid_steps = ['all', '1', '2', '3', '4', '5a', '5b', '6', '7']
s_flag = 0
for step in steps:
    if step not in valid_steps:
        sys.stderr.write('Invalid step detected: ' + step + '\n')
        s_flag = 1
if s_flag:
    sys.stderr.write('Invalid step or steps used.  Please specify all or list as csv, i.e 1,2,3 (no spaces!)\n')
    sys.exit(1)

config_file = args.config_file
init_in = [config_file]

if config_file is not None:
    sys.stderr.write('Doing a quick check to see if config file ' + config_file + ' exists\n')
else:
    sys.stderr.write('No config file set!  Need file for param -j or --config\n')
    sys.exit(1)
check_files(init_in, 'WRAP', '0')
with open(config_file) as f:
    config_data = json.load(f)
script_dir = config_data['script_dir']
if script_dir[-1] != '/':
    script_dir += '/'
sys.stderr.write('Doing a quick check to see if script dir ' + script_dir + ' exists\n')
check_dirs([script_dir], 'WRAP', '0')
sys.stderr.write('Doing a quick check to see if relevant vars in config file are set and files and dirs actually '
                 'exist for step(s) ' + args.steps + '\n')
config_precheck()

if 'all' in steps or '1' in steps:
    step1 = '1_query_cavatica.py'
    cmd = script_dir + step1 + ' --output ' + args.cav + ' --profile ' + args.profile + ' --config ' \
          + config_file + ' 2> step1.errs'
    sys.stderr.write('Running step 1: ' + cmd + '\n')
    sys.stderr.flush()
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as e:
        print(e)
        sys.stderr.write('Failed at step 1.  Maybe check cavatica profile and if sevenbridges api in installed?\n')
        sys.exit(1)
    outputs = [args.cav]
    check_files(outputs, 'OUTPUT', '1')

if 'all' in steps or '2' in steps:
    step2 = '2_query_ds_by_bs_id.py'
    inputs = [args.cav]
    check_files(inputs, 'INPUT', step2)
    cmd = script_dir + step2 + ' --kf-url ' + args.url + ' --cavatica ' + args.cav + ' --config ' \
          + config_file + ' 2> step2.errs'
    sys.stderr.write('Running step 2: ' + cmd + '\n')
    sys.stderr.flush()
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as e:
        print(e)
        sys.stderr.write('Failed at step 2. Maybe the url is wrong or inaccessible?\n')
        sys.exit(1)
    outputs = ['tum_bs_ds_info.txt', 'norm_bs_ds_info.txt']
    check_files(outputs, 'OUTPUT', '2')

if 'all' in steps or '3' in steps:
    tumor = args.tumor
    normal = args.normal
    step3 = '3_create_data_sheets.py'
    if '2' in steps or 'all' in steps:
        tumor = 'tum_bs_ds_info.txt'
        normal = 'norm_bs_ds_info.txt'
    inputs = [tumor, normal, args.cav]
    check_files(inputs, 'INPUT', step3)
    cmd = script_dir + step3 + ' --cavatica ' + args.cav + ' --config ' + config_file + ' --tumor ' \
          + tumor + ' --normal ' + normal
    if args.cl_supp is not None:
        cmd += ' --cell-line-supplement ' + args.cl_supp
    cmd += ' 2> step3.errs'
    sys.stderr.write('Running step 3: ' + cmd + '\n')
    sys.stderr.flush()
    subprocess.call(cmd, shell=True)
    outputs = ['datasheets']
    check_dirs(outputs, 'OUTPUT', '3')

if 'all' in steps or '4' in steps:
    step4 = '4_merge_maf.py'
    input_files = [args.header, args.cav]
    check_files(input_files, 'INPUT', step4)
    input_dirs = [args.maf_dir]
    check_dirs(input_dirs, 'INPUT', step4)
    cmd = script_dir + step4 + ' --header ' + args.header + ' --cavatica ' + args.cav + ' --config ' \
          + config_file + ' --maf-dir ' + args.maf_dir + ' 2> step4.errs'
    sys.stderr.write('Running step 4: ' + cmd + '\n')
    sys.stderr.flush()
    subprocess.call(cmd, shell=True)
    output_dirs = ['merged_mafs']
    check_dirs(output_dirs, 'OUTPUT', '4')

if 'all' in steps or '5a' in steps:
    step5a = '5a_cnv_genome2gene.py'
    if config_data['cna_flag']:
        input_dirs = [args.cnv_dir]
        check_dirs(input_dirs, 'INPUT', step5a)
        cmd = script_dir + step5a + ' --config ' + config_file + ' --cnv-dir ' + args.cnv_dir \
              + ' 2> step5a.errs'
        sys.stderr.write('Running step 5a: ' + cmd + '\n')
        sys.stderr.flush()
        try:
            subprocess.check_call(cmd, shell=True)
        except Exception as e:
            print(e)
            sys.stderr.write('Failed at step 5a.  Maybe bedtools is improperly defined?\n')
            sys.exit(1)
        output_dirs = ['converted_cnvs']
        check_dirs(output_dirs, 'OUTPUT', '5a')
    else:
        sys.stderr.write('cna_flag set to 0. Check that you have CNA data to process and set to 1 in the config file\n')
        sys.stderr.flush()

if 'all' in steps or '5b' in steps:
    step5b = '5b_merge_cnv.py'
    if config_data['cna_flag']:
        input_files = [args.cav]
        check_files(input_files, 'INPUT', step5b)
        cnv_gene_dir = args.cnv_dir_gene
        if 'all' in steps or '5a' in steps:
            cnv_gene_dir = 'converted_cnvs'
        input_dirs = [cnv_gene_dir]
        check_dirs(input_dirs, 'INPUT', step5b)
        cmd = script_dir + step5b + '  --cavatica ' + args.cav + ' --config ' + config_file + ' --cnv-dir-gene ' \
              + cnv_gene_dir + ' 2> step5b.errs'
        sys.stderr.write('Running step 5b: ' + cmd + '\n')
        sys.stderr.flush()
        subprocess.call(cmd, shell=True)
        output_dirs = ['merged_cnvs']
        check_dirs(output_dirs, 'OUTPUT', '5b')
    else:
        sys.stderr.write('cna_flag set to 0. Check that you have CNA data to process and set to 1 in the config file\n')
        sys.stderr.flush()

if 'all' in steps or '6' in steps:
    step6 = '6_merge_rsem.py'
    if config_data['rna_flag']:
        input_files = [args.cav]
        check_files(input_files, 'INPUT', step6)
        input_dirs = [args.rsem_dir]
        check_dirs(input_dirs, 'INPUT', step6)
        cmd = script_dir + step6 + '  --cavatica ' + args.cav + ' --config ' + config_file + ' --rsem-dir ' \
              + args.rsem_dir + ' --profile ' + args.profile + ' 2> step6.errs'
        sys.stderr.write('Running step 6: ' + cmd + '\n')
        sys.stderr.flush()
        subprocess.call(cmd, shell=True)
        output_dirs = ['merged_rsem']
        check_dirs(output_dirs, 'OUTPUT', '6')
    else:
        sys.stderr.write('rna_flag set to 0. Check that you have RNA data to process and set to 1 in the config file\n')
        sys.stderr.flush()

if 'all' in steps or '7' in steps:
    input_files = []
    if config_data['rna_flag']:
        input_files.append(config_data['fpkmFileLoc'])
        check_files(input_files, 'INPUT', '7')
    cmd = 'Rscript ' + script_dir + '7_createLoadingFiles.R ' + config_file + ' 2> step7.errs > step7.out'
    sys.stderr.write('Running step 7: ' + cmd + '\n')
    sys.stderr.flush()
    subprocess.call(cmd, shell=True)
    output_dirs = ['processed']
    check_dirs(output_dirs, 'OUTPUT', '7')

sys.stderr.write('Pipe complete, check logs\n')
