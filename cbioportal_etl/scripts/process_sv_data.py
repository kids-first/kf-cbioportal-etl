#!/usr/bin/env python3
"""Converts our current fusion and annot SV inputs into the cBio SV format.

For repeat rows (same breakpoint in a sample, different callers),
ARRIBA annot is used with ceiling of mean for counts used
"""

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd
from annotsv_to_cbio_sv import init_cbio_sv_df, sv_setup_outdir_metadata
from convert_fusion_as_sv import fusion_setup_outdir_metadata, init_cbio_master
import pdb


def main():
    parser = argparse.ArgumentParser(
        description="Convert openPBTA fusion table OR list of annofuse files to cbio format."
    )
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        help="Table with cbio project, kf bs ids, cbio IDs, and file names",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--fusion-results",
        action="store",
        help="annoFuse results dir OR openX merged fusion file",
        required=False,
    )
    parser.add_argument(
        "--sv-results",
        action="store_true",
        help="DNA SV results from annotSV",
        required=False,
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        action="store",
        dest="out_dir",
        default="merged_sv",
        help="Result output dir. Default is merged_fusion",
    )
    parser.add_argument(
        "-m",
        "--mode",
        action="store",
        dest="mode",
        help="describe source, openX or kfprod or dgd",
        required=False,
    )
    parser.add_argument(
        "-a",
        "--append",
        action=argparse.BooleanOptionalAction,
        dest="append",
        help="Flag to append, meaning print to STDOUT and skip header",
        required=False,
    )

    args = parser.parse_args()
    if args.fusion_results is None and args.sv_results is None:
        print(
            "Either --fusion-results and/or --sv-results must be specified.",
            file=sys.stderr,
        )
        sys.exit(1)

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    sv_results, fusion_results = args.sv_results, args.fusion_results
    if sv_results:
        print("DNA SV flag given, processing", file=sys.stderr)
        # get DNA SV results if any
        link_input_dir = "annotSV_results"
        dna_sv_subset: pd.DataFrame = sv_setup_outdir_metadata(link_input_dir, args.out_dir, args.table)
        cbio_sv_df: pd.DataFrame = init_cbio_sv_df(link_input_dir, dna_sv_subset)
    if fusion_results:
        print("RNA fusion data given, processing", file=sys.stderr)
        # Get RNA Fusion results if any
        if args.mode is None:
            print("Need to set -m mode when processing fusion data. Opts are kfprod, openX, dgd", file=sys.stderr)
            sys.exit(1)
        # ext used in pbta vs openpedcan varies
        rna_fusion_subset: pd.DataFrame = fusion_setup_outdir_metadata(args.mode, args.out_dir, args.table)

        cbio_fusion_df = init_cbio_master(fusion_results, args.mode, rna_fusion_subset)

    if sv_results and fusion_results:
        print("Both DNA SV and RNA fusion given, merging data frames", file=sys.stderr)
        merged_metadata_subset = pd.concat([dna_sv_subset, rna_fusion_subset])
        project_list: np.ndarray = merged_metadata_subset.cbio_project.unique()
        merged_df: pd.DataFrame = pd.concat([cbio_sv_df, cbio_fusion_df])
    elif sv_results:
        merged_metadata_subset = dna_sv_subset
        project_list: np.ndarray = dna_sv_subset.cbio_project.unique()
        merged_df = cbio_sv_df
    else:
        merged_metadata_subset = rna_fusion_subset
        project_list: np.ndarray = dna_sv_subset.cbio_project.unique()
        merged_df = cbio_sv_df

    for project in project_list:
        print(f"Outputting final SV results for {project}", file=sys.stderr)
        sub_sample_list = list(
            merged_metadata_subset.loc[merged_metadata_subset["cbio_project"] == project, "cbio_sample_name"]
        )
        out_fname = os.path.join(args.out_dir, project + ".svs.txt")
        subset_merged_df: pd.DataFrame = merged_df[merged_df.Sample_Id.isin(sub_sample_list)]
        #subset_merged_df.fillna("", inplace=True)
        if not args.append:
            subset_merged_df.set_index("Sample_Id", inplace=True)
            subset_merged_df.to_csv(out_fname, sep="\t", mode="w", index=True, quoting=csv.QUOTE_NONE)
        else:
            # append to existing fusion file, use same header
            existing = pd.read_csv(out_fname, sep="\t", keep_default_na=False, na_values=[""])
            subset_merged_df = subset_merged_df[existing.columns]
            subset_merged_df.set_index("Sample_Id", inplace=True)
            subset_merged_df.to_csv(
                out_fname, sep="\t", mode="a", index=True, quoting=csv.QUOTE_NONE, header=None
            )
if __name__ == "__main__":
    main()
