#!/usr/bin/env python3

import sys
import argparse
import os
import subprocess
import pdb

def process_col_samps(h_cnt, s_col, fname, samp_dict):
    cur = open(fname)
    temp = {}
    for x in range(h_cnt):
        skip = next(cur)
    for row in cur:
        data = row.rstrip('\n').split('\t')
        samp = data[s_col]
        if samp not in temp:
            temp[samp] = 1
            if samp not in samp_dict:
                samp_dict[samp] = 1
            else:
                samp_dict[samp] += 1
    cur.close()
    return samp_dict

def process_row_samps(s_col, fname, samp_dict):
    cur = open(fname)
    head = next(cur)
    header = head.rstrip('\n').split('\t')
    for x in range(s_col, len(header), 1):
        if header[x] not in samp_dict:
            samp_dict[header[x]] = 1
        else:
            samp_dict[header[x]] += 1
    cur.close()
    return samp_dict


parser = argparse.ArgumentParser(description='Add case lists to existing projects')
parser.add_argument('-t', '--type', action='store', dest='ctype', help='case_list type to gen, i.e. cnaseq, 3way_complete, SCRIPT SHOULD BE RUN RIGHT AT STUDY DIR LEVEL')

args = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

# files to search for case building
clist_fdict = {'cnaseq': {'fnames': ['data_mutations_extended.txt', 'data_CNA.txt'], 
            'name': 'Tumor samples with mutatation and CNA data','desc': 'All tumor samples with mutation and CNA data','cat': 'all_cases_with_mutation_and_cna_data'} , 
            '3way_complete': { 'fnames': ['data_mutations_extended.txt', 'data_CNA.txt','data_rna_seq_v2_mrna.txt'], 
            'name': 'Tumor samples with mutatation, CNA and mRNA data','desc': 'All tumor samples with mutation, CNA, and mRNA data','cat': 'all_cases_with_mutation_and_cna_and_mrna_data'}}

f_cmd = 'find . -name ' + clist_fdict[args.ctype]['fnames'][0]
result = subprocess.check_output(f_cmd, shell=True)
muts_list = result.decode().split('\n')
for mutFile in muts_list:
    if mutFile != '':
        s_dict = {}
        project = os.path.dirname(mutFile)
        project = project.replace('./','')
        sys.stderr.write("Generating case lists for project " + project + "\n")
        sys.stderr.write("Processing " + mutFile + "\n")
        s_dict = process_col_samps(2, 13, mutFile, s_dict)
        for i in range(1, len(clist_fdict[args.ctype]['fnames'])):
            fname = project + "/" +  clist_fdict[args.ctype]['fnames'][i]
            sys.stderr.write("Processing " + fname + "\n")
            s_dict = process_row_samps(1, fname, s_dict)

        new_case_list = open(project + "/case_lists/cases_" + args.ctype + ".txt", 'w')
        new_case_list.write('cancer_study_identifier: ' + project + '\nstable_id: ' + project + '_' + args.ctype + '\ncase_list_name: ' + clist_fdict[args.ctype]['name'] + '\n')
        count = len(clist_fdict[args.ctype]['fnames'])
        s_list = []
        for samp in s_dict:
            if s_dict[samp] == count:
                s_list.append(samp)
        total = len(s_list)
        new_case_list.write('case_list_description: ' + clist_fdict[args.ctype]['desc'] + ' (' + str(total) + ' samples)\ncase_list_category: ' + clist_fdict[args.ctype]['cat'] + '\ncase_list_ids: ' + '\t'.join(s_list))
        new_case_list.close()
        sys.stderr.write("Completed creating case list for " + project + "\n")
    else:
        sys.stderr.write(clist_fdict[args.ctype]['fnames'][0] + " path was empty str.  skpping!\n")
