#!/usr/bin/env python3

"""
Script to append dgd fusion results to pbta using STDOUT or collate if standalone
"""
import sys
import os
import argparse
from get_file_metadata_helper import get_file_metadata

parser = argparse.ArgumentParser(description='Output fields DGD fusion - meant to be appended to an existing file!')
parser.add_argument('-f', '--fusion-dir', action='store', dest='fusion_dir', help='Fusion file directory')
parser.add_argument('-t', '--table', action='store', dest='table',
                    help='Table with cbio project, kf bs ids, cbio IDs, and file names')
parser.add_argument('-a', '--append', action='store_true', dest='append',
                    help='Optional - if given will output to stdout to append, else will create new merged file')
parser.add_argument('-o', '--out-dir', action='store', dest='out_dir', default='merged_fusion/', help='Result output dir. Default is merged_fusion')



args = parser.parse_args()
file_meta_dict = get_file_metadata(args.table, 'DGD_FUSION')

header = ['Hugo_Symbol','Entrez_Gene_Id','Center','Tumor_Sample_Barcode','Fusion','DNA_support','RNA_support','Method','Frame']
if not args.append:
    try:
        os.mkdir(args.out_dir)
    except:
        sys.stderr.write('output dir already exists\n')
        sys.stderr.flush()
mid_idx = header.index('Method')
t_idx = header.index('Tumor_Sample_Barcode')
for cbio_dx in file_meta_dict:
    if not args.append:
        out_file = open(args.out_dir + cbio_dx + '.fusions.txt')
    else:
        out_file = sys.stdout
    for cbio_tum_id in file_meta_dict[cbio_dx]:
        fusion = file_meta_dict[cbio_dx][cbio_tum_id]['fname']
        sys.stderr.write("Processing " + fusion + "\n")
        cur = open(args.fusion_dir + "/" + fusion)
        f_head = next(cur)
        for data in cur:
            datum = data.rstrip('\n').split('\t')
            # Set method
            datum[mid_idx] = 'DGD_curated'
            datum[t_idx] = cbio_tum_id
            out_file.write("\t".join(datum) + "\n")
        sys.stderr.write("Processed " + fusion + "\n")
    out_file.close()