#!/usr/bin/env python3

import sys
import argparse
import os
import pdb

parser = argparse.ArgumentParser(description='Added survival data to data_clinical_patient.txt sheet. Run from dir level of disease subdirs')
parser.add_argument('-d', '--data', action='store', dest='data', help='File with patients IDS and survival data')
parser.add_argument('-l', '--list', action='store', dest='flist', help='Data file list')
parser.add_argument('-o', '--outdir', action='store', dest='outdir', help='output dir')
parser.add_argument('-t', '--id-type', action='store', dest='itype', help='kf or ext')
parser.add_argument('-f', '--dfs-flag', action='store', dest='dflag', help='1 if DFS data available, 0 if not')

args = parser.parse_args()

added_header = ['OS_STATUS\tOS_MONTHS\tDFS_STATUS\tDFS_MONTHS',
    'Overall patient survival status\tOverall survival in months since initial diagnosis\tDisease free status since initial treatment\tDisease free (months) since initial treatment',
    'STRING\tNUMBER\tSTRING\tNUMBER',
    '1\t1\t1\t1',
    'OS_STATUS\tOS_MONTHS\tDFS_STATUS\tDFS_MONTHS']
if int(args.dflag) == 0:
    added_header = ['OS_STATUS\tOS_MONTHS',
    'Overall patient survival status\tOverall survival in months since initial diagnosis',
    'STRING\tNUMBER',
    '1\t1',
    'OS_STATUS\tOS_MONTHS']

if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)
surv_dict = {}
surv_data = open(args.data)
skip_head = next(surv_data)
for line in surv_data:
    data = line.rstrip('\n').split('\t')
    surv_dict[data[0]] = '\t'.join(data[1:])
surv_data.close()
id_idx = {'kf': 0, 'ext': 1}
for fname in open(args.flist):
    fname = fname.rstrip('\n')
    subdir = os.path.dirname(fname)
    newdir = args.outdir + '/' + subdir
    if not os.path.isdir(newdir):
        os.makedirs(newdir, exist_ok=True)
    cur = open(fname)
    newdata = open(newdir + '/data_clinical_patient.txt', 'w')
    # pdb.set_trace()
    for i in range(0, len(added_header)):
        head = next(cur)
        newdata.write(head.rstrip('\n') + '\t' + added_header[i] + '\n')
    for line in cur:
        newdata.write(line.rstrip('\n'))
        data = line.split('\t')
        pt_id = data[id_idx[args.itype]]
        if pt_id in surv_dict:
            newdata.write('\t' + surv_dict[pt_id] + '\n')
        else:
            newdata.write('\t\t\t\t\n')
            sys.stderr.write('No survival entry for ' + pt_id + ' in ' + fname + '\n')
    cur.close()
    newdata.close()
