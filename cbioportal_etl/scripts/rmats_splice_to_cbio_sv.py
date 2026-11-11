#!/usr/bin/env python3
"""Converts rMATs inputs into the cBio SV format.

Multiple files per sample comprise different splice types
"""

import argparse
import csv
import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
import pandas as pd
from annotsv_to_cbio_sv import sv_setup_outdir_metadata

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)



def eval_effect_on_frame(row: pd.Series) -> str:
    """Evaluate event info for each row of rmats file.

    Args:
        row: Pandas series with rmats row data

    Returns:
        Event info string for cBio SV file

    """
    site2_effect_on_frame = "NA"
    if row["Class"] != "RETAINED_INTRON":
        site2_effect_on_frame = "In_frame" if not row["SV_Length"] % 3  else "Frameshift"
    return site2_effect_on_frame



def rename_and_format_cols(concat_sv_df: pd.DataFrame) -> pd.DataFrame:
    """Rename and format columns for cBio SV file.

    Args:
        concat_sv_df: Pandas dataframe with concatenated rmats data
    Returns:
        Formatted dataframe for cBio SV file
    """
    logger = logging.getLogger(__name__)
    logger.info("Renaming and formatting annotSV columns for cBio SV output")
    desired: list[str] = [
        "Sample_Id",
        "SV_Status",
        "Site1_Hugo_Symbol",
        "Site1_Chromosome",
        "Site1_Position",
        "Site1_Region",
        "Site2_Hugo_Symbol",
        "Site2_Chromosome",
        "Site2_Position",
        "Site2_Region",
        "Site2_Effect_On_Frame",
        "NCBI_Build",
        "Class",
        "Tumor_Read_Count",
        "Event_Info",
        "Annotation",
        "DNA_Support",
        "SV_Length",
        "Comments",
    ]

    concat_sv_df["Site2_Effect_On_Frame"] = concat_sv_df.apply(eval_effect_on_frame, axis=1)
    # need to add custom RMATS. cols for IJC, SJC, IncLevel, copy flanking base pair info to Annotation


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

    # track two values - one is part of file name, other is desired column value for cBio SV file
    event_types: dict[str, list[str]] = {"alternative_3_prime_splice_sites": ["filtered.A3SS", "ALT_3SS"],
                         "alternative_5_prime_splice_sites": ["filtered.A5SS", "ALT_5SS"],
                         "retained_introns": ["filtered.RI", "RETAINED_INTRON"],
                         "skipped_exons": ["filtered.SE", "SKIPPED_EXON"],
                         }
    frame_list = []
    for metadata in splice_metadata.itertuples():
        # get event type from file_type, as columns slightly different and will inform cbio entries
        event_type = None
        for _file_type, event_code in event_types.items():
            if event_code[0] in metadata.file_name:
                event_type = event_code[1]
                break
        if event_type is None:
            value_err = f"Could not find valid event type for {metadata.file_name}"
            raise ValueError(value_err)
        # concat rmats files, rename Sample Column according to cBio name
        try:
            fname = f"{splice_results}/{metadata.file_name}"
            rmats_file: pd.DataFrame = pd.read_csv(
                fname,
                sep="\t", keep_default_na=False, na_values=[""],
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
        tumor_read_count = rmats_file["IJC_SAMPLE_1"] + rmats_file["SJC_SAMPLE_1"]
        sv_length = site2_position - site1_position

        # assignment fest
        rmats_file = rmats_file.assign(Sample_Id=metadata.cbio_sample_name,
                                       SV_Status="SOMATIC",
                                       Site1_Hugo_Symbol = rmats_file["geneSymbol"],
                                       Site1_Chromosome=rmats_file["chr"],
                                       Site1_Position=site1_position,
                                       Site1_Region=region,
                                       Site2_Hugo_Symbol = rmats_file["geneSymbol"],
                                       Site2_Chromosome=rmats_file["chr"],
                                       Site2_Position=site2_position,
                                       Site2_Region=region,
                                       Class=event_type,
                                       Tumor_Read_Count=tumor_read_count,
                                       SV_Length=sv_length,
                                       RNA_Support="Yes",
                                       NCBI_Build="GRCh38"
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
    cbio_sv_df = rename_and_format_cols(concat_sv_df=cbio_sv_df)

if __name__ == "__main__":
    main()
