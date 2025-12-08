#!/usr/bin/env python3
"""Merge MAF Files.

Merge and filter MAF files to create an output with a singular header and
variants not nor,ally suitable for cBio removed

"""

import argparse
import concurrent.futures
import json
import os
import sys
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
    maf_fn: str, new_maf: IO, maf_exc: dict[str, int], tum_id: str, norm_id: str
) -> None:
    """Iterate over maf file, skipping header lines since the files are being merged.

    With possibility of mixed source, search headers

    Args:
    maf_fn: Input MAF filename
    new_maf: Output maf file handle
    maf_exc: Dict with Variant_Classification exclusionary terms
    tum_id: ID to change Tumor_Sample_Barcode to
    norm_id: ID to change Matched_Norm_Sample_Barcode to

    """
    with open(maf_fn) as cur_maf:
        next(cur_maf)
        head: str = next(cur_maf)

        cur_header: list = head.rstrip("\n").split("\t")
        h_dict: dict[str, int | None] = {}
        for item in print_header:
            if item in cur_header:
                h_dict[item] = cur_header.index(item)
            else:
                h_dict[item] = None

        entry_indices = EntryIndices(
            cur_header.index("Tumor_Sample_Barcode"),
            cur_header.index("Matched_Norm_Sample_Barcode"),
            cur_header.index("Variant_Classification"),
            cur_header.index("Hugo_Symbol"),
        )

        with concurrent.futures.ThreadPoolExecutor(16) as executor:
            results: set[Future[list[str] | None]] = {
                executor.submit(
                    filter_entry,
                    entry,
                    tum_id,
                    norm_id,
                    entry_indices,
                    maf_exc,
                )
                for entry in cur_maf
            }
            for result in concurrent.futures.as_completed(results):
                filtered: list[str] | None = result.result()
                if filtered is not None:
                    to_print: list[str] = []
                    for item in print_header:
                        if h_dict[item] is not None:
                            to_print.append(filtered[h_dict[item]])
                        else:
                            to_print.append("")

                    new_maf.write("\t".join(to_print) + "\n")


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
                process_maf(maf_dir + fname, new_maf, maf_exc, cbio_tum_id, cbio_norm_id)
                x += 1
            print(f"Completed processing {x} entries in {study}", file=sys.stderr)
            new_maf.close()
    except Exception as e:
        print(e, file=sys.stderr)
        sys.exit()


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
