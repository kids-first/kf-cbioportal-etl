#!/usr/bin/env python3
import sys
import os
import argparse
import json
from get_file_metadata_helper import get_file_metadata
import concurrent.futures


def filter_entry(entry, tum_id, norm_id, tid_idx, nid_idx, v_idx, eid_idx):
    """
    Only output entries not in exclusion list while dropping ENTREZ ID, but keeping TERT promoter hits
    """
    data = entry.rstrip("\n").split("\t")
    # Want to allow TERT promoter as exception to exclusion rules
    if data[v_idx] not in maf_exc or (
        data[h_idx] == "TERT" and data[v_idx] == "5'Flank"
    ):
        data[tid_idx] = tum_id
        data[nid_idx] = norm_id
        data.pop(eid_idx)
        return data
    else:
        return None


def process_maf(
    maf_fn, new_maf, maf_exc, tum_id, norm_id, tid_idx, h_idx, nid_idx, v_idx, eid_idx
):
    """
    Iterate over maf file, skipping header lines since the files are being merged
    """
    cur_maf = open(maf_fn)
    next(cur_maf)
    next(cur_maf)
    with concurrent.futures.ThreadPoolExecutor(16) as executor:
        results = {
            executor.submit(
                filter_entry, entry, tum_id, norm_id, tid_idx, nid_idx, v_idx, eid_idx
            )
            for entry in cur_maf
        }
        for result in concurrent.futures.as_completed(results):
            if result.result() != None:
                new_maf.write("\t".join(result.result()) + "\n")
    cur_maf.close()


def process_tbl(
    cbio_dx, file_meta_dict, tid_idx, h_idx, nid_idx, v_idx, eid_idx, print_head
):
    """
    Probaby a less likely scenario, but can split out into multiple projects based on dict
    """
    try:
        x = 0
        # project/disease name should be name of directory hosting datasheet
        sys.stderr.write("Processing " + cbio_dx + " project" + "\n")
        new_maf = open(out_dir + cbio_dx + ".maf", "w")
        new_maf.write(print_head)
        for cbio_tum_id in file_meta_dict[cbio_dx]:
            cbio_norm_id = file_meta_dict[cbio_dx][cbio_tum_id]["cbio_norm_id"]
            fname = file_meta_dict[cbio_dx][cbio_tum_id]["fname"]
            sys.stderr.write(
                "Found relevant maf to process for "
                + " "
                + cbio_tum_id
                + " "
                + cbio_norm_id
                + " "
                + file_meta_dict[cbio_dx][cbio_tum_id]["kf_tum_id"]
                + " "
                + file_meta_dict[cbio_dx][cbio_tum_id]["kf_norm_id"]
                + " "
                + fname
                + "\n"
            )
            sys.stderr.flush()
            process_maf(
                maf_dir + fname,
                new_maf,
                maf_exc,
                cbio_tum_id,
                cbio_norm_id,
                tid_idx,
                h_idx,
                nid_idx,
                v_idx,
                eid_idx,
            )
            x += 1
        sys.stderr.write(
            "Completed processing " + str(x) + " entries in " + cbio_dx + "\n"
        )
        new_maf.close()
    except Exception as e:
        print(e)
        sys.exit()


parser = argparse.ArgumentParser(
    description="Merge and filter mafs using cavatica task info and datasheets."
)
parser.add_argument(
    "-t",
    "--table",
    action="store",
    dest="table",
    help="Table with cbio project, kf bs ids, cbio IDs, and file names",
)
parser.add_argument(
    "-i", "--header", action="store", dest="header", help="File with maf header only"
)
parser.add_argument(
    "-m", "--maf-dir", action="store", dest="maf_dir", help="maf file directory"
)
parser.add_argument(
    "-j",
    "--config",
    action="store",
    dest="config_file",
    help="json config file with data types and " "data locations",
)
parser.add_argument(
    "-f",
    "--dgd-status",
    action="store",
    dest="dgd_status",
    help="Flag to determine load will have pbta/kf + dgd(both), kf/pbta only(kf), dgd-only(dgd)",
    default="both",
    const="both",
    nargs="?",
    choices=["both", "kf", "dgd"],
)


args = parser.parse_args()
with open(args.config_file) as f:
    config_data = json.load(f)
# get maf file ext
maf_dir = args.maf_dir
if maf_dir[-1] != "/":
    maf_dir += "/"
# If DGD maf only, else if both, dgd maf wil be handled separately, or no at all if no dgd and kf only
file_meta_dict = get_file_metadata(args.table, "DGD_MAF")
if args.dgd_status != "dgd":
    file_meta_dict = get_file_metadata(args.table, "maf")
head_fh = open(args.header)

print_head = next(head_fh)
cur_head = next(head_fh)
cur_header = cur_head.rstrip("\n").split("\t")
eid_idx = cur_header.index("Entrez_Gene_Id")
tid_idx = cur_header.index("Tumor_Sample_Barcode")
nid_idx = cur_header.index("Matched_Norm_Sample_Barcode")
v_idx = cur_header.index("Variant_Classification")
h_idx = cur_header.index("Hugo_Symbol")
cur_header.pop(eid_idx)

head_fh.close()
print_head += "\t".join(cur_header) + "\n"
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
out_dir = "merged_mafs/"
try:
    os.mkdir(out_dir)
except:
    sys.stderr.write("output dir already exists\n")

for cbio_dx in file_meta_dict:
    process_tbl(
        cbio_dx, file_meta_dict, tid_idx, h_idx, nid_idx, v_idx, eid_idx, print_head
    )

sys.stderr.write("Done, check logs\n")
