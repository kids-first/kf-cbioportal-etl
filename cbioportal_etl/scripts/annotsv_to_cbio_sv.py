#!/usr/bin/env python3
"""Converts our current DNA annotSV inputs into the cBio SV format.

Focuses on the "split" rows of the Annotation_mode to get per-gene SV data
"""

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd
import pdb


def setup_outdir_metadata(link_input_dir: str, out_dir: str, table: str) -> pd.DataFrame:
    """Create output dir and link input dir if they don't exist, then subset metadata for DNA SV files.

    Output dir will store the merged SV files, and link_input_dir will store symlinks to the annotSV results for processing.

    Args:
        link_input_dir: Directory to store symlinks to annotSV results
        out_dir: Directory to store merged SV results
        table: Table with cbio project, kf bs ids, cbio IDs, and file names
    Returns:
        Pandas dataframe with DNA SV-specific metadata

    """
    # create output dir and link input dir if they don't exist
    os.makedirs(link_input_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    # deal only with SV metadata
    ext = "sv"
    # ensure sample name is imported as str
    all_file_meta: pd.DataFrame = pd.read_csv(table, sep="\t", dtype={"cbio_sample_name": str})
    dna_sv_subset: pd.DataFrame = all_file_meta.loc[all_file_meta["etl_file_type"] == ext]
    # Create symlinks to mafs in one place for ease of processing
    sv_dirs_in: str = ",".join(dna_sv_subset.file_type.unique().tolist())
    print(f"Symlinking SV files from {sv_dirs_in} to {link_input_dir}", file=sys.stderr)
    sym_errs: int = 0
    for dirname in sv_dirs_in.split(","):
        if os.path.exists(dirname):
            abs_path: str = os.path.abspath(dirname)
        else:
            print(f"Input SV dir {dirname} does not exist. Check if files were downloaded to the correct location", file=sys.stderr)
            sys.exit(2)
        try:
            for fname in os.listdir(dirname):
                src: str = os.path.join(abs_path, fname)
                dest: str = os.path.join(link_input_dir, fname)
                os.symlink(src, dest)
        except Exception as e:
            print(e, file=sys.stderr)
            print(f"Could not sym link {fname} in {dirname}", file=sys.stderr)
            sym_errs += 1
    # If symlink errors, stop here as data will be incomplete
    if sym_errs:
        print(f"Could not sym link {sym_errs} files, exiting!", file=sys.stderr)
        sys.exit(2)
    return dna_sv_subset


def init_cbio_sv_df(sv_results: str, sv_metadata: pd.DataFrame) -> pd.DataFrame:
    """Use data frame subset on DNA SV files to find and merge result files.

    Args:
        sv_results: annotSV results dir
        sv_metadata: Pandas dataframe with DNA SV-specific metadata
    Returns:
        Collapsed and formatted dataframe, list of desired fields

    """
    desired: list[str] = [
        "SV_chrom",
        "SV_start",
        "SV_length",
        "SV_type",
        "INFO",
        "FORMAT",
        "CytoBand",
        "Gene_name",
        "Frameshift",
    ]
    flist: pd.Series[str] = sv_metadata.file_name
    frame_list = []
    try:
        for i in range(0, len(flist), 1):
            # concat annot SV file, rename Sample Column according to cBio name
            ann_file: pd.DataFrame = pd.read_csv(
                f"{sv_results}/{flist.iloc[i]}", sep="\t", keep_default_na=False, na_values=[""]
            )
            # drop entries where Annotation_mode is not "split" (i.e. per-gene)
            ann_file: pd.DataFrame = ann_file.loc[ann_file["Annotation_mode"] == "split"]
            ann_file = ann_file.assign(Sample=sv_metadata.iloc[i].cbio_sample_name)
            # at this step, get the tumor BS ID
            # parse FORMAT for PR and SR and assign values cBio columns
            # drop the tumor and normal BS ID cols
            affected_id: str = sv_metadata.iloc[i].affected_bs_id
            reference_id: str = sv_metadata.iloc[i].reference_bs_id
            read_support = ann_file[affected_id].str.split(":", expand=True)
            ann_file["Tumor_Paired_End_Read_Count"] = (
                read_support[0].str.split(",").str[1].astype("Int64")
            )
            # some files have 0 SR entries at all, so we need to handle that case
            if 1 in read_support.columns:
                ann_file["Tumor_Split_Read_Count"] = pd.to_numeric(
                read_support[1].str.split(",").str[1],
                errors="coerce",
                ).astype("Int64")
            else:
                ann_file["Tumor_Split_Read_Count"] = pd.Series(pd.NA, index=ann_file.index, dtype="Int64")
            ann_file = ann_file.drop(columns=[affected_id, reference_id])
            frame_list.append(ann_file)
    except Exception as e:
        print(f"{e}", file=sys.stderr)
        pdb.set_trace()
        sys.exit(1)
    concat_frame: pd.DataFrame = pd.concat(frame_list)
    pdb.set_trace()
    hold = 1
    return concat_frame
    # concat_frame = filter_and_format_annots(sample_renamed_df=concat_frame, drop_low=True)
    # del frame_list
    # fusion_data: pd.DataFrame = concat_frame[desired]


def main():
    parser = argparse.ArgumentParser(
        description="Convert annotSV results to cbio format."
    )
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        help="Table with cbio project, kf bs ids, cbio IDs, and file names",
        required=True,
    )
    parser.add_argument(
        "-l",
        "--link-input-dir",
        action="store",
        help="DNA SV directory name to store annotSV symlinks",
        default="annotSV_results/",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        action="store",
        dest="out_dir",
        default="merged_sv/",
        help="Result output dir. Default is merged_sv",
    )
    args = parser.parse_args()

    dna_sv_subset: pd.DataFrame = setup_outdir_metadata(args.link_input_dir, args.out_dir, args.table)
    pdb.set_trace()
    cbio_sv_df: pd.DataFrame = init_cbio_sv_df(args.link_input_dir, dna_sv_subset)


if __name__ == "__main__":
    main()