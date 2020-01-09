#!/usr/bin/env python3

import sys
import argparse
import os

parser = argparse.ArgumentParser(description='Swaps data between two projects. Run from dir level of disease subdirs')
parser.add_argument('-a', '--datainput1', action='store', dest='d1', help='First sheet with sample to swap')
parser.add_argument('-b', '--datainput2', action='store', dest='d2', help='Second sheet with sample to swap')
parser.add_argument('-x', '--sample1', action='store', dest='x', help='sample 1 name to swap with...')
parser.add_argument('-y', '--sample2', action='store', dest='y', help='sample2 name to swap')
parser.add_argument('-o', '--outdir', action='store', dest='outdir', help='outputdir')
parser.add_argument('-f', '--flag', action='store', dest='f', help='0 to keep original sampe names, 1 to swap '
                                                                      'names too')

args = parser.parse_args()
d1_out_dir = args.outdir + '/' + os.path.dirname(args.d1)
d2_out_dir = args.outdir + '/' + os.path.dirname(args.d2)
if not os.path.isdir(d1_out_dir):
    os.mkdir(d1_out_dir)
if not os.path.isdir(d2_out_dir):
    os.mkdir(d2_out_dir)

d1_out_fh = open(d1_out_dir + '/' + os.path.basename(args.d1), 'w')
d2_out_fh = open(d2_out_dir + '/' + os.path.basename(args.d2), 'w')

d1_in = open(args.d1)
d2_in = open(args.d2)
# find sample one to swap
h1 = next(d1_in)
h1_info = h1.rstrip('\n').split('\t')
h1_i = 0
for i in range(0, len(h1_info)):
    if h1_info[i] == args.x:
        h1_i = i
        sys.stderr.write('Found sample ' + args.x + ' to swap in ' + args.d1 + ' at position ' + str(i) + '\n')
        if args.f == '1':
            h1_info[i] = args.y
            sys.stderr.write('Swapped header name to ' + args.y + '\n')
        break
d1_out_fh.write('\t'.join(h1_info) + '\n')

#find sample 2 to swap
h2 = next(d2_in)
h2_info = h2.rstrip('\n').split('\t')
h2_i = 0
for i in range(0, len(h2_info)):
    if h2_info[i] == args.y:
        h2_i = i
        sys.stderr.write('Found sample ' + args.y + ' to swap in ' + args.d2 + ' at position ' + str(i) + '\n')
        if args.f == '1':
            h2_info[i] = args.x
            sys.stderr.write('Swapped header name to ' + args.x + '\n')
        break
d2_out_fh.write('\t'.join(h2_info) + '\n')

for line1 in d1_in:
    line2 = next(d2_in)

    info1 = line1.rstrip('\n').split('\t')
    info2 = line2.rstrip('\n').split('\t')
    t1 = info1[h1_i]
    info1[h1_i] = info2[h2_i]
    info2[h2_i] = t1
    d1_out_fh.write('\t'.join(info1) + '\n')
    d2_out_fh.write('\t'.join(info2) + '\n')
d1_out_fh.close()
d2_out_fh.close()