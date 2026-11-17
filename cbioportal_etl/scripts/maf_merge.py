#!/usr/bin/env python3
"""Merge MAF Files.

Merge and filter MAF files to create an output with a singular header and
variants not normally suitable for cBio removed

"""

import argparse
import json
import os
import sys
from operator import itemgetter
from typing import IO, NamedTuple

import pandas as pd


class EntryIndices(NamedTuple):
    """Named tuple to simplify passing on indicies of headers."""

    tid_idx: int
    nid_idx: int
    v_idx: int
    h_idx: int


def filter_entry(
    entry: str,
    tum_id: str,
    norm_id: str,
    entry_indices: EntryIndices,
    maf_exc: dict[str, int],
) -> list[str] | None:
    """Only output entries not in exclusion list while dropping ENTREZ ID, but keeping TERT promoter hits.

    Args:
    entry: line from input MAF file
    tum_id: ID to change Tumor_Sample_Barcode to
    norm_id: ID to change Matched_Norm_Sample_Barcode to
    entry_indices: Named tuple with teh following:
        Index position for Tumor_Sample_Barcode ID in entry as list,
        Index position for Matched_Norm_Sample_Barcode ID in entry as list,
        Index position of Variant_Classification in entry as list,
        Index position of Hugo_Symbol in entry as list
    maf_exc: Dictionary of variant classifications to exclude from output

    Returns:
    list of entry values with Tumor_Sample_Barcode and Matched_Norm_Sample_Barcode changed to
      cbio IDs and ENTREZ ID dropped, or None if entry is to be excluded

    """
    data: list = entry.rstrip("\n").split("\t")
    # Want to allow TERT promoter as exception to exclusion rules
    if data[entry_indices.v_idx] not in maf_exc or (
        data[entry_indices.h_idx] == "TERT" and data[entry_indices.v_idx] == "5'Flank"
    ):
        data[entry_indices.tid_idx] = tum_id
        data[entry_indices.nid_idx] = norm_id
        return data
    return None


def process_maf(
        maf_fn: str,
        new_maf: IO,
        maf_exc: dict[str, int],
        tum_id: str,
        norm_id: str,
        print_header: list[str]) -> None:
    """Process a single MAF file, filtering out unwanted entries and writing to new MAF."""
    # read in ~ 1MB at a time
    with open(maf_fn, buffering=1 << 20) as cur_maf:
        next(cur_maf)  # first line, skipped as in the original
        cur_header = next(cur_maf).rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(cur_header)}
        tid_idx = col["Tumor_Sample_Barcode"]
        nid_idx = col["Matched_Norm_Sample_Barcode"]
        v_idx = col["Variant_Classification"]
        h_idx = col["Hugo_Symbol"]

        # Output columns missing from this file get index -1, which points at the
        # "" appended to each row below.
        pick = itemgetter(*[col.get(name, -1) for name in print_header])
        split_n = max(v_idx, h_idx) + 1
        write = new_maf.write

        for line in cur_maf:
            # Cheap partial split first: only the leading fields are needed to decide.
            head = line.split("\t", split_n)
            vc = head[v_idx]
            if vc in maf_exc and not (vc == "5'Flank" and head[h_idx] == "TERT"):
                continue
            row = line.rstrip("\n").split("\t")
            row[tid_idx] = tum_id
            row[nid_idx] = norm_id
            row.append("")
            write("\t".join(pick(row)) + "\n")


def process_tbl(study: str, file_meta_dict: list[list], print_head: str) -> None:
    """Process MAF file be project (cbio_dx).

    Unlikely that more than one project will be processed, pass each variant entry for filtering
    Args:
    cbio_dx: cBio project name
    file_meta_dict: Dict that has been subset by file type from ETL file
    print_head: Output header line for file output an dto guide output content
    """
    try:
        x = 0
        # project/disease name should be name of directory hosting datasheet
        print(f"Processing {study} project", file=sys.stderr)
        with open(f"{out_dir}{study}.maf", "w") as new_maf:
            new_maf.write(print_head)
            for fname, cbio_tum_id, cbio_norm_id in file_meta_dict:
                print(
                    f"Found relevant maf to process for {cbio_tum_id} {cbio_norm_id} {fname}",
                    file=sys.stderr,
                )
                process_maf(maf_dir + fname, new_maf, maf_exc, cbio_tum_id, cbio_norm_id, print_header)
                x += 1
            print(f"Completed processing {x} entries in {study}", file=sys.stderr)
    except Exception as e:
        print(e, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge and filter MAFs using ETL table withe file locations.",
    )
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        dest="table",
        help="Table with cbio project, kf bs ids, cbio IDs, and file names",
    )
    parser.add_argument(
        "-i",
        "--header",
        action="store",
        dest="header",
        help="File with maf header only",
    )
    parser.add_argument(
        "-j",
        "--config",
        action="store",
        dest="config_file",
        help="json config file with data types and data locations",
    )

    args = parser.parse_args()
    with open(args.config_file) as f:
        config_data = json.load(f)
    # Create symlinks to mafs in one place for ease of processing
    maf_dir: str = "MAFS/"
    all_file_meta = pd.read_csv(args.table, sep="\t", dtype={"cbio_sample_name": str})
    maf_subset = all_file_meta[all_file_meta["etl_file_type"] == "maf"].copy()
    study = maf_subset["cbio_project"].iloc[0]
    maf_list = maf_subset[["file_name", "cbio_sample_name", "cbio_matched_normal_name"]].fillna("").astype(str).drop_duplicates().values.tolist()

    maf_dirs_in: str = ",".join(maf_subset.file_type.unique().tolist())
    print(f"Symlinking maf files from {maf_dirs_in} to {maf_dir}", file=sys.stderr)
    os.makedirs("MAFS", exist_ok=True)
    sym_errs: int = 0
    for dirname in maf_dirs_in.split(","):
        if os.path.exists(dirname):
            abs_path: str = os.path.abspath(dirname)
        else:
            print(f"Input MAF dir {dirname} does not exist. Check if files were downloaded to the correct location", file=sys.stderr)
            sys.exit(2)
        try:
            for fname in os.listdir(dirname):
                src: str = os.path.join(abs_path, fname)
                dest: str = os.path.join(maf_dir, fname)
                os.symlink(src, dest)
        except Exception as e:
            print(e, file=sys.stderr)
            print(f"Could not sym link {fname} in {dirname}", file=sys.stderr)
            sym_errs += 1
    # If symlink errors, stop here as data will be incomplete
    if sym_errs:
        print(f"Could not sym link {sym_errs} files, exiting!", file=sys.stderr)
        sys.exit(2)

    with open(args.header) as head_fh:
        print_head: str = next(head_fh)
        cur_head: str = next(head_fh)
        print_header: list[str] = cur_head.rstrip("\n").split("\t")
        eid_idx: int = print_header.index("Entrez_Gene_Id")
        print_header.pop(eid_idx)

    print_head += "\t".join(print_header) + "\n"
    maf_exc: dict[str, int] = {
        "Silent": 0,
        "Intron": 0,
        "IGR": 0,
        "3'UTR": 0,
        "5'UTR": 0,
        "3'Flank": 0,
        "5'Flank": 0,
        "RNA": 0,
    }
    out_dir: str = "merged_mafs/"
    os.makedirs(out_dir, exist_ok=True)
    # iterating through projects that are the first key in dict
    process_tbl(study, maf_list, print_head)

    sys.stderr.write("Done, check logs\n")
