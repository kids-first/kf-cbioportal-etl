#!/usr/bin/env python3

"""
Script to append dgd fusion results to pbta using STDOUT or collate if standalone
"""
import sys
import os
import argparse
import gzip
import pdb
from get_file_metadata_helper import get_file_metadata


parser = argparse.ArgumentParser(
    description="Output fields DGD fusion - meant to be appended to an existing file!"
)
parser.add_argument(
    "-f",
    "--fusion-dir",
    action="store",
    dest="fusion_dir",
    help="Fusion file directory",
)
parser.add_argument(
    "-t",
    "--table",
    action="store",
    dest="table",
    help="Table with cbio project, kf bs ids, cbio IDs, and file names",
)
parser.add_argument(
    "-a",
    "--append",
    action="store_false",
    dest="append",
    help="Optional - if given will output to stdout to append, else will create new merged file",
)
parser.add_argument(
    "-m",
    "--merged",
    action="store_true",
    dest="merged",
    help="If input is already merged, treat fusion dir as file instead",
)
parser.add_argument(
    "-o",
    "--out-dir",
    action="store",
    dest="out_dir",
    default="merged_fusion/",
    help="Result output dir. Default is merged_fusion",
)

args = parser.parse_args()

header = [
    "Hugo_Symbol",
    "Entrez_Gene_Id",
    "Center",
    "Tumor_Sample_Barcode",
    "Fusion",
    "DNA_support",
    "RNA_support",
    "Method",
    "Frame",
]
if args.append:
    try:
        os.mkdir(args.out_dir)
    except:
        sys.stderr.write("output dir already exists\n")
        sys.stderr.flush()
mid_idx = header.index("Method")
t_idx = header.index("Tumor_Sample_Barcode")
if args.merged:
    # merged file needs to get cBio ID by BS ID instead of using table to parse individual files and assin cBio IDs
    file_meta_dict = {}
    meta_table =  open(args.table)
    mhead = next(meta_table)
    mheader = mhead.rstrip('\n').split('\t')
    mb_idx = mheader.index('T_CL_BS_ID')
    cbio_idx = mheader.index('Cbio_Tumor_Name')
    for line in meta_table:
        info = line.rstrip('\n').split('\t')
        file_meta_dict[info[mb_idx]] = info[cbio_idx]
else:
    file_meta_dict = get_file_metadata(args.table, "DGD_FUSION")

if args.merged:
    out_file = sys.stdout
    if args.fusion_dir[-3:] == '.gz':
        cur = gzip.open(args.fusion_dir, mode='rt')
    else:
        cur = open(args.fusion_dir)
    f_head = next(cur)
    for data in cur:
        datum = data.rstrip("\n").split("\t")
        # To fit current format of having a blank value for Entrez ID column
        datum.insert(1,"")
        cbio_tum_id = file_meta_dict[datum[t_idx]]
        # Set method
        datum[mid_idx] = "DGD_curated"
        datum[t_idx] = cbio_tum_id
        out_file.write("\t".join(datum) + "\n")
else:
    for cbio_dx in file_meta_dict:
        if args.append:
            out_file = open(args.out_dir + cbio_dx + ".fusions.txt", "w")
            out_file.write("\t".join(header) + "\n")
        else:
            out_file = sys.stdout
        for cbio_tum_id in file_meta_dict[cbio_dx]:
            fusion = file_meta_dict[cbio_dx][cbio_tum_id]["fname"]
            sys.stderr.write("Processing " + fusion + "\n")
            cur = open(args.fusion_dir + "/" + fusion)
            f_head = next(cur)
            for data in cur:
                datum = data.rstrip("\n").split("\t")
                # Set method
                datum[mid_idx] = "DGD_curated"
                datum[t_idx] = cbio_tum_id
                out_file.write("\t".join(datum) + "\n")
            sys.stderr.write("Processed " + fusion + "\n")
        out_file.close()
