#!/usr/bin/env python3
"""
Helper script to append DGD data to an existing merged maf file.
Uses filter_entry to filter out undesired calls like in other mafs
"""
import sys
import argparse
from get_file_metadata_helper import get_file_metadata
from maf_merge import filter_entry

parser = argparse.ArgumentParser(
    description="Output fields from maf file based on header - meant to be appended to an existing file!"
)
parser.add_argument(
    "-i", "--header", action="store", dest="header", help="File with maf header only"
)
parser.add_argument(
    "-m", "--maf-dir", action="store", dest="maf_dir", help="maf file directory"
)
parser.add_argument(
    "-t",
    "--table",
    action="store",
    dest="table",
    help="Table with cbio project, kf bs ids, cbio IDs, and file names",
)

args = parser.parse_args()
file_meta_dict = get_file_metadata(args.table, "DGD_MAF")

header_file = open(args.header)
maf_version = next(header_file)
head = next(header_file)
# store header fields as list, will use to get position in new file, if exists
# will output blank if not in new file
header = head.rstrip("\n").split("\t")
eid_idx = header.index("Entrez_Gene_Id")
header.pop(eid_idx)
h_dict = {}
# dict of classifications to drop
maf_exc = {
    "Silent": 0,
    "Intron": 0,
    "IGR": 0,
    "3'UTR": 0,
    "5'UTR": 0,
    "3'Flank": 0,
    "5'Flank": 0,
    "RNA": 0,
}

for item in header:
    h_dict[item] = None
# Set filler for norm ID
norm_id=""
for cbio_dx in file_meta_dict:
    skipped = 0
    for cbio_tum_id in file_meta_dict[cbio_dx]:
        maf = file_meta_dict[cbio_dx][cbio_tum_id]["fname"]
        sys.stderr.write("Processing " + maf + "\n")
        cur = open(args.maf_dir + "/" + maf)
        maf_version = next(cur)
        m_head = next(cur)
        m_header = m_head.rstrip("\n").split("\t")
        tid_idx = m_header.index("Tumor_Sample_Barcode")
        nid_idx = m_header.index("Matched_Norm_Sample_Barcode")
        v_idx = m_header.index("Variant_Classification")
        h_idx = m_header.index("Hugo_Symbol")

        for i in range(len(m_header)):
            if m_header[i] in h_dict:
                h_dict[m_header[i]] = i
        # only print items in original header and in same order, else print blank
        for data in cur:
            to_print = []
            datum = filter_entry(data, cbio_tum_id, norm_id, tid_idx, nid_idx, v_idx, h_idx, maf_exc)
            # Set tumor barcode to cBio ID
            if datum:
                for item in header:
                    if h_dict[item] != None:
                        to_print.append(datum[h_dict[item]])
                    else:
                        to_print.append("")
                sys.stdout.write("\t".join(to_print) + "\n")
            else:
                skipped += 1
        sys.stderr.write("Processed " + maf + "\n")
        sys.stderr.write("Skipped " + str(skipped) + " entries meeting exlusion criteria\n")
