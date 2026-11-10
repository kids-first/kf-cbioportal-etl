#!/usr/bin/env python3
"""Converts rMATs inputs into the cBio SV format.

Multiple files per sample comprise different splice types
"""

import argparse
import csv
from importlib.metadata import metadata
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
    logger.info("Processing %d rmats files from %s", len(splice_metadata), splice_results)

    event_types: dict = {"alternative_3_prime_splice_sites": "ALT_3SS",
                         "alternative_5_prime_splice_sites": "ALT_5SS",
                         "retained_introns": "RETAINED_INTRON",
                         "skipped_exons": "SKIPPED_EXON",
                         }
    frame_list = []
    for metadata in splice_metadata.itertuples():
        # get event type from file_type, as columns slightly different and will inform cbio entries
        event_type = None
        for file_type, event_code in event_types.items():
            if file_type in metadata.file_name:
                event_type = event_code
                break
        if event_type is None:
            value_err = f"Could not find valid event type for {metadata.file_name}"
            raise ValueError(value_err)
        # concat rmats files, rename Sample Column according to cBio name
        try:
            fname = f"{splice_results}/{metadata.file_name}"
            rmats_file: pd.DataFrame = pd.read_csv(
                fname,
                sep="\t", keep_default_na=False, na_values=[""]
            )
        except Exception as e:
            file_process_err = "Could not open files %s in %s", metadata.file_name, splice_results
            raise FileNotFoundError(file_process_err) from e
        # set event type-specific values
        region = "INTRON" if event_type == "RETAINED_INTRON" else "EXON"
        site1_position = rmats_file["longExonStart_0base"] if event_type in ("ALT_3SS", "ALT_5SS") \
            else rmats_file["exonStart_0base"] if event_type == "SKIPPED_EXON" \
                else rmats_file["riExonStart_0base"]
        site2_position = rmats_file["longExonEnd"] if event_type in ("ALT_3SS", "ALT_5SS") \
                    else rmats_file["exonEnd"] if event_type == "SKIPPED_EXON" \
                        else rmats_file["riExonEnd"]
        sv_length = site2_position.to_numpy() - site1_position.to_numpy()
        site2_effect_on_frame = "NA"
        if event_type != "RETAINED_INTRON":
            site2_effect_on_frame = "In_frame" if not sv_length % 3  else "Frameshift"
        # assignment fest
        rmats_file = rmats_file.assign(Sample_Id=metadata.cbio_sample_name,
                                       Site1_Hugo_Symbol = rmats_file["geneSymbol"],
                                       Site1_Chromosome=rmats_file["chr"],
                                       Site1_Position=site1_position,
                                       Site1_Region=region,
                                       Site2_Hugo_Symbol = rmats_file["geneSymbol"],
                                       Site2_Chromosome=rmats_file["chr"],
                                       Site2_Position=site2_position,
                                       Site2_Region=region,
                                       Class=event_type,
                                       Site2_Effect_on_Frame=site2_effect_on_frame,
                                       )

        frame_list.append(rmats_file)
    concat_frame: pd.DataFrame = pd.concat(frame_list)
    del frame_list
    return concat_frame


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
