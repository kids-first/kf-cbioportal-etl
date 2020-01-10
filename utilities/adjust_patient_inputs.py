#!/usr/bin/env python3

import sys
import argparse
import os
import re
from shutil import copyfile
import concurrent.futures
# import pdb


def process_case_files(old_cdir, new_cdir, case_files, rm_sample_list):
    f = 0
    for fn in case_files:
        old_cf = open(old_cdir + '/' + fn)
        new_cf = open(new_cdir + '/' + fn, 'w')
        new_s_list = []
        # print standard lines
        standard = next(old_cf)
        new_cf.write(standard)
        standard = next(old_cf)
        new_cf.write(standard)
        standard = next(old_cf)
        new_cf.write(standard)
        # hold lines to change description
        desc = next(old_cf)
        buffer = next(old_cf)
        entry = next(old_cf)
        info = entry.rstrip('\n').split()
        for i in range(1, len(info)):
            if info[i] not in rm_sample_list:
                new_s_list.append(info[i])
        f = len(new_s_list)
        desc = re.sub('\d+', str(f), desc)
        new_cf.write(desc + buffer + info[0] + ' ' + '\t'.join(new_s_list) + '\n')
        old_cf.close()
        new_cf.close()

    return f


def process_data_clinical(old_ds, new_ds, rm_patient_list):
    for entry in old_ds:
        info = entry.split('\t')
        if info[0] not in rm_patient_list:
            new_ds.write(entry)


def process_samp_in_head(old_file, new_file, rm_sample_list):
    try:
        old_head = next(old_file)
        good_i = []
        head_split = old_head.rstrip('\n').split('\t')
        new_head = []
        new_head.append(head_split[0])
        for i in range(1, len(head_split)):
            if head_split[i] not in rm_sample_list:
                good_i.append(i)
                new_head.append(head_split[i])
        new_file.write('\t'.join(new_head) + '\n')
        for data in old_file:
            datum = data.rstrip('\n').split('\t')
            new_out = datum[0]
            for i in good_i:
                new_out += '\t' + datum[i]
            new_file.write(new_out + '\n')
    except Exception as e:
        print(e)
        exit(1)


def process_samp_in_tbl(old_file, new_file, rm_sample_list):
    # pdb.set_trace()
    head = next(old_file)
    new_file.write(head)
    for data in old_file:
        datum = data.split('\t')
        if datum[13] not in rm_sample_list:
            new_file.write(data)


def cleanup_all(dname, rm_sample_list, rm_patient_list, new_out):
    try:
        # redo case list files
        new_cdir = new_out + '/case_lists'
        os.mkdir(new_cdir)
        old_cdir = dname + '/case_lists'
        case_files = os.listdir(old_cdir)
        n_samp = process_case_files(old_cdir, new_cdir, case_files, rm_sample_list)
        sys.stderr.write('Reprocessed case lists for ' + dname + '\n')
        sys.stderr.flush()
        if n_samp == 0:
            sys.stderr.write('All samples in ' + dname + ' were REMOVED!\n')
        else:
            # redo data_clinical sheets
            new_dss = open(new_out + '/data_clinical_sample.txt', 'w')
            old_dss = open(dname + '/data_clinical_sample.txt')
            process_data_clinical(old_dss, new_dss, rm_patient_list)
            new_dss.close()

            new_dps = open(new_out + '/data_clinical_patient.txt', 'w')
            old_dps = open(dname + '/data_clinical_patient.txt')
            process_data_clinical(old_dps, new_dps, rm_patient_list)
            new_dps.close()
            sys.stderr.write('Reprocessed data clinical sheets for ' + dname + '\n')
            sys.stderr.flush()

            samp_in_head_fn = ['data_CNA.txt', 'data_linear_CNA.txt', 'data_rna_seq_v2_mrna_median_Zscores.txt',
                               'data_rna_seq_v2_mrna.txt']
            for fn in samp_in_head_fn:
                old_fh = open(dname + '/' + fn)
                new_fh = open(new_out + '/' + fn, 'w')
                process_samp_in_head(old_fh, new_fh, rm_sample_list)
                new_fh.close()
                sys.stderr.write('Reprocessed ' + fn + ' for ' + dname + '\n')
                sys.stderr.flush()

            samp_in_tbl = ['data_mutations_extended.txt']
            for fn in samp_in_tbl:
                old_fh = open(dname + '/' + fn)
                new_fh = open(new_out + '/' + fn, 'w')
                process_samp_in_tbl(old_fh, new_fh, rm_sample_list)
                new_fh.close()
                sys.stderr.write('Reprocessed ' + fn + ' for ' + dname + '\n')
                sys.stderr.flush()
            # copy over all meta tables
            meta_list = os.listdir(dname)
            for fn in meta_list:
                if fn[0:5] == 'meta_':
                    copyfile((dname + '/' + fn), (new_out + '/' + fn))
                    sys.stderr.write('Copied meta file ' + fn + '\n')
                    sys.stderr.flush()
    except Exception as e:
        print(e)


def check_data_sheet(dname, exc_list, outdir):
    try:
        cur = open(dname + '/data_clinical_sample.txt')
        f = 0
        rm_slist = []
        rm_plist = []
        for line in cur:
            info = line.rstrip('\n').split('\t')
            if info[0] in exc_list:
                f = 1
                rm_slist.append(info[1])
                if info[0] not in rm_plist:
                    rm_plist.append(info[0])
        cur.close()
        if f == 0:
            sys.stderr.write('No samples to remove from ' + dname + '\n')
            sys.stderr.flush()
        else:
            sys.stderr.write(dname + ' has samples ' + str(len(rm_slist)) + ' to be excluded:' + '\t'.join(rm_slist)
                             + ' from patients ' + '\t'.join(rm_plist) + '\n')
            sys.stderr.flush()
            parts = dname.split('/')
            new_out = outdir + '/' + parts[-1]
            os.mkdir(new_out)
            cleanup_all(dname, rm_slist, rm_plist, new_out)
        return dname
    except Exception as e:
        print(e)
        exit(1)


parser = argparse.ArgumentParser(description='Use list of sample IDs to remove to rebuild a project.')
parser.add_argument('-l', '--list', action='store', dest='plist', help='new-line separated patient list')
parser.add_argument('-d', '--dir', action='store', dest='dir', help='data dir to search')
parser.add_argument('-p', '--pattern', action='store', dest='pattern', help='dir search pattern - a suffix or prefix')
parser.add_argument('-o', '--outdir', action='store', dest='outdir', help='output directory for projects that change')

args = parser.parse_args()

# read in patient list
exc_list = []

for line in open(args.plist):
    exc_list.append(line.rstrip('\n'))

# get dir list
all_files = os.listdir(args.dir)
dir_list = []
pattern = '.*' + args.pattern + '.*'
for fn in all_files:
    if re.match(pattern, fn):
        dir_list.append(args.dir + '/' + fn)
# for pdir in dir_list:
#     check_data_sheet(pdir, exc_list, args.outdir)
with concurrent.futures.ProcessPoolExecutor(8) as executor:
    results = {executor.submit(check_data_sheet, pdir, exc_list, args.outdir): pdir for pdir in dir_list}
    for result in concurrent.futures.as_completed(results):
        sys.stderr.write('Completed processing ' + result.result() + '\n')
