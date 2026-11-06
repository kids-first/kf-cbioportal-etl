#!/usr/bin/env python3
"""Converts rMATs inputs into the cBio SV format.

Multiple files per sample comprise different splice types
"""

import argparse
import csv
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from annotsv_to_cbio_sv import sv_setup_outdir_metadata

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)



def init_cbio_sv_df(splice_results: str, splice_metadata: pd.DataFrame) -> pd.DataFrame:
    """Use data frame subset on DNA SV files to find and merge result files.

    Args:
        splice_results: rMATs splice results dir
        splice_metadata: Pandas dataframe with rMATs splice-specific metadata
    Returns:
        Collapsed and formatted dataframe, list of desired fields

    """
    logger = logging.getLogger(__name__)
    logger.info("Processing %d annotSV files from %s", len(splice_metadata), splice_results)
    frame_list = []
    for metadata in splice_metadata.itertuples():
        # concat annot SV file, rename Sample Column according to cBio name
        try:
            fname = f"{splice_results}/{metadata.file_name}"
            ann_file: pd.DataFrame = pd.read_csv(
                fname,
                sep="\t", keep_default_na=False, na_values=[""]
            )
        except Exception as e:
            file_process_err = "Could not open files %s in %s", metadata.file_name, splice_results
            raise FileNotFoundError(file_process_err) from e
    return cbio_sv_df

def main() -> None:
    """Convert rmats splice results to cBio SV format."""
    parser = argparse.ArgumentParser(
        description="Convert rmats splice results to cbio format.",
    )
    parser.add_argument(
        "-t",
        "--table",
        help="Table with cbio project, kf bs ids, cbio IDs, and file names",
        required=True,
    )
    parser.add_argument(
        "-l",
        "--link-input-dir",
        help="rmats splice directory name to store symlinks",
        default="rmats_results/",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        default="merged_sv/",
        help="Result output dir. Default is merged_sv",
    )
    parser.add_argument(
        "-s",
        "--skip-symlink",
        action="store_true",
        help="Flag to skip sym links if repeating the process with the same input files",
        required=False,
    )

    args = parser.parse_args()

    rna_splice_subset: pd.DataFrame = sv_setup_outdir_metadata(args.link_input_dir, args.out_dir,
                                                            args.table, "splice", args.skip_symlink)
    project_list: np.ndarray = rna_splice_subset.cbio_project.unique()
    cbio_sv_df: pd.DataFrame = init_cbio_sv_df(args.link_input_dir, rna_splice_subset)
