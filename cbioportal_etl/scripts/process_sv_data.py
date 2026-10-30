#!/usr/bin/env python3
"""Converts our current fusion and annot SV inputs into the cBio SV format.

For repeat rows (same breakpoint in a sample, different callers),
ARRIBA annot is used with ceiling of mean for counts used
"""

import argparse
import csv
import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
import pandas as pd
from annotsv_to_cbio_sv import init_cbio_sv_df, rename_and_format_cols, sv_setup_outdir_metadata
from convert_fusion_as_sv import (
    fusion_setup_outdir_metadata,
    init_cbio_master,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def main() -> None:
    """Convert fusion and/or annotSV results to cBio SV format."""
    logger = logging.getLogger(__name__)
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
    )
    parser.add_argument(
        "--process-sv-results",
        action=argparse.BooleanOptionalAction,
        help="DNA SV results from annotSV exist in table and you wish to process them",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        action="store",
        default="merged_sv",
        help="Result output dir. Default is merged_fusion",
    )
    parser.add_argument(
        "-m",
        "--mode",
        action="store",
        choices=["openX", "kfprod", "dgd"],
        help="describe source",
    )
    parser.add_argument(
        "-a",
        "--append",
        action=argparse.BooleanOptionalAction,
        help="Flag to append, meaning print to STDOUT and skip header",
    )
    parser.add_argument(
        "-s",
        "--skip-symlink",
        action=argparse.BooleanOptionalAction,
        help="Flag to skip sym links if repeating the process with the same input files",
    )

    args = parser.parse_args()
    # ensure dependent args are set
    if not (args.fusion_results or args.process_sv_results):
        parser.error("Either --fusion-results and/or --sv-results must be specified.")
    if args.fusion_results and not args.mode:
        parser.error("Need to set -m mode when processing fusion data. Opts are kfprod, openX, dgd.")

    out_dir = args.out_dir
    Path(out_dir).mkdir(exist_ok=True, parents=True)

    if args.process_sv_results:
        logger.info("DNA SV flag given, processing")
        # get DNA SV results if any
        link_input_dir = "annotSV_results"
        dna_sv_subset: pd.DataFrame = sv_setup_outdir_metadata(link_input_dir, args.out_dir, args.table, args.skip_symlink)
        cbio_sv_df: pd.DataFrame = init_cbio_sv_df(link_input_dir, dna_sv_subset)
        cbio_sv_df = rename_and_format_cols(concat_sv_df=cbio_sv_df)
        merged_metadata_subset = dna_sv_subset
        project_list: np.ndarray = dna_sv_subset.cbio_project.unique()
        merged_df = cbio_sv_df

    if args.fusion_results:
        logger.info("RNA fusion data given, processing")
        # Get RNA Fusion results if any
        rna_fusion_subset: pd.DataFrame = fusion_setup_outdir_metadata(args.mode, args.out_dir, args.table)

        cbio_fusion_df = init_cbio_master(args.fusion_results, args.mode, rna_fusion_subset)
        merged_metadata_subset = rna_fusion_subset
        merged_df = cbio_fusion_df

    if args.process_sv_results and args.fusion_results:
        logger.info("Both DNA SV and RNA fusion given, merging data frames")
        merged_metadata_subset = pd.concat([dna_sv_subset, rna_fusion_subset])
        merged_df: pd.DataFrame = pd.concat([cbio_sv_df, cbio_fusion_df])

    project_list: np.ndarray = merged_metadata_subset.cbio_project.unique()
    for project in project_list:
        logger.info("Outputting final SV results for %s", project)
        sub_sample_list = list(
            merged_metadata_subset.loc[merged_metadata_subset["cbio_project"] == project, "cbio_sample_name"]
        )
        out_fname = Path(args.out_dir) / (project + ".svs.txt")
        subset_merged_df: pd.DataFrame = merged_df[merged_df.Sample_Id.isin(sub_sample_list)]
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
