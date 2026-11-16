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



EVENT_TYPES = {  # file name token -> cBio Class
    "filtered.A3SS": "ALT_3SS",
    "filtered.A5SS": "ALT_5SS",
    "filtered.RI": "RETAINED_INTRON",
    "filtered.SE": "SKIPPED_EXON",
}
ANNOTATION_COLS = {
    "ALT_3SS": ["shortES", "shortEE", "flankingES", "flankingEE"],
    "ALT_5SS": ["shortES", "shortEE", "flankingES", "flankingEE"],
    "RETAINED_INTRON": ["upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "SKIPPED_EXON": ["upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
}
POSITION_COLS = {  # (site1 start, site2 end)
    "ALT_3SS": ("longExonStart_0base", "longExonEnd"),
    "ALT_5SS": ("longExonStart_0base", "longExonEnd"),
    "SKIPPED_EXON": ("exonStart_0base", "exonEnd"),
    "RETAINED_INTRON": ("riExonStart_0base", "riExonEnd"),
}
BASE_COLS = ["geneSymbol", "chr", "strand", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1"]


def make_annotation_event_info(df: pd.DataFrame, event_type: str) -> tuple[pd.Series, pd.Series]:
    """Build Annotation and Event_Info series (no mutation of df)."""
    annotation = pd.Series("", index=df.index, dtype=object)
    for col in ANNOTATION_COLS[event_type]:
        annotation += f"{col}: " + df[col].astype(str) + ";"
    event_info = (
        df["Site1_Hugo_Symbol"].astype(str) + f"_{event_type}_"
        + df["chr"].astype(str) + ":"
        + df["Site1_Position"].astype(str) + "-" + df["Site2_Position"].astype(str)
    )
    return annotation, event_info


def populate_annotation_event_info(rmats_file_df: pd.DataFrame,
                                   event_type: str) -> tuple[pd.Series, pd.Series, list[str]]:
    """Populate annotation column based on event type.

    Args:
        rmats_file_df: Pandas dataframe with rmats file data
        event_type: Type of splice event

    Returns:
        Series with populated annotation column, series with populated event info column, and list of cols to drop later

    """
    cols_for_annotate: dict[str, list[str]] = {
        "ALT_3SS": ["shortES", "shortEE", "flankingES", "flankingEE"],
        "ALT_5SS": ["shortES", "shortEE", "flankingES", "flankingEE"],
        "RETAINED_INTRON": ["upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
        "SKIPPED_EXON": ["upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    }
    annotation = pd.Series("", index=rmats_file_df.index, dtype=object)
    for col in cols_for_annotate[event_type]:
        annotation += f"{col}: " + rmats_file_df[col].astype(str) + ";"
    event_info = (
        rmats_file_df["Site1_Hugo_Symbol"].astype(str)
        + "_" + rmats_file_df["Class"].astype(str)
        + "_" + rmats_file_df["chr"].astype(str)
        + ":" + rmats_file_df["Site1_Position"].astype(str)
        + "-" + rmats_file_df["Site2_Position"].astype(str)
    )
    return annotation, event_info, cols_for_annotate[event_type]


def rename_and_format_cols(concat_sv_df: pd.DataFrame, comment: str) -> pd.DataFrame:
    """Rename and format columns for cBio SV file.

    Args:
        concat_sv_df: Pandas dataframe with concatenated rmats data
        comment: Comment to add to the comments column

    Returns:
        Formatted dataframe for cBio SV file

    """
    logger = logging.getLogger(__name__)
    logger.info("Renaming and formatting rMATS dataframe columns for cBio SV output")
    desired: list[str] = [
        "cbio_project",
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
        "RNA_Support",
        "SV_Length",
        "Comments",
        "RMATS.strand",
        "RMATS.IJC_SAMPLE_1",
        "RMATS.SJC_SAMPLE_1",
        "RMATS.IncLevel1",
    ]
    rename_dict: dict[str, str] = {
        "strand": "RMATS.strand",
        "IJC_SAMPLE_1": "RMATS.IJC_SAMPLE_1",
        "SJC_SAMPLE_1": "RMATS.SJC_SAMPLE_1",
        "IncLevel1": "RMATS.IncLevel1",
    }

    df = concat_sv_df.rename(columns=rename_dict)
    df["Site2_Effect_On_Frame"] = np.where(
        df["Class"] == "RETAINED_INTRON",
        "NA",
        np.where(df["SV_Length"] % 3 == 0, "In_frame", "Frameshift"),
    )
    return df.assign(Comments=comment)[desired]


def init_cbio_sv_df(splice_results: str, splice_metadata: pd.DataFrame) -> pd.DataFrame:
    """Use data frame subset on rMATS splice files to find and merge result files.

    Args:
        splice_results: rMATs splice results dir
        splice_metadata: Pandas dataframe with rMATs splice-specific metadata
    Returns:
        Collapsed and formatted dataframe, list of desired fields

    """
    logger = logging.getLogger(__name__)
    logger.info("Processing %d rmats files from %s", len(splice_metadata), splice_results)

    frame_list = []
    n = 25
    for i, metadata in enumerate(splice_metadata.itertuples(), start=1):
        event_type = next(
            (code for token, code in EVENT_TYPES.items() if token in metadata.file_name), None
        )
        if event_type is None:
            raise ValueError(f"Could not find valid event type for {metadata.file_name}")

        start_col, end_col = POSITION_COLS[event_type]
        usecols = [*BASE_COLS, start_col, end_col, *ANNOTATION_COLS[event_type]]
        fname = Path(splice_results) / metadata.file_name
        try:
            rmats_file = pd.read_csv(
                fname, sep="\t", keep_default_na=False, na_values=[""], usecols=usecols
            )
        except (OSError, ValueError, pd.errors.ParserError) as e:
            raise RuntimeError(f"Could not read {fname}: {e}") from e

        site1_position = rmats_file[start_col]
        site2_position = rmats_file[end_col]
        rmats_file = rmats_file.assign(
            cbio_project=metadata.cbio_project,
            Sample_Id=metadata.cbio_sample_name,
            SV_Status="SOMATIC",
            Site1_Hugo_Symbol=rmats_file["geneSymbol"],
            Site1_Chromosome=rmats_file["chr"],
            Site1_Position=site1_position,
            Site1_Region="INTRON" if event_type == "RETAINED_INTRON" else "EXON",
            Site2_Hugo_Symbol=rmats_file["geneSymbol"],
            Site2_Chromosome=rmats_file["chr"],
            Site2_Position=site2_position,
            Site2_Region="INTRON" if event_type == "RETAINED_INTRON" else "EXON",
            Class=event_type,
            Tumor_Read_Count=rmats_file["IJC_SAMPLE_1"] + rmats_file["SJC_SAMPLE_1"],
            SV_Length=site2_position - site1_position,
            RNA_Support="Yes",
            NCBI_Build="GRCh38",
        )
        rmats_file["Annotation"], rmats_file["Event_Info"] = make_annotation_event_info(
            rmats_file, event_type
        )
        rmats_file = rmats_file.drop(
            columns=[start_col, end_col, "geneSymbol", "chr", *ANNOTATION_COLS[event_type]]
        )
        frame_list.append(rmats_file)
        if i % n == 0:
            logger.info("Processed %d of %d rmats files", i, len(splice_metadata))
    concat_frame: pd.DataFrame = pd.concat(frame_list, ignore_index=True)
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
    parser.add_argument(
        "--comment",
        help="rMATS and GENCODE version information to add to the comments column in the cBio SV file",
        default="rMATS v4.1.2, GENCODE v39",
    )

    args = parser.parse_args()

    rna_splice_subset: pd.DataFrame = sv_setup_outdir_metadata(args.link_input_dir, args.out_dir,
                                                            args.table, "splice", args.skip_symlink)
    cbio_sv_df: pd.DataFrame = init_cbio_sv_df(args.link_input_dir, rna_splice_subset)
    comment = args.comment
    cbio_sv_df = rename_and_format_cols(cbio_sv_df, comment)

    for project, proj_df in cbio_sv_df.groupby("cbio_project"):
        out = Path(args.out_dir) / f"{project}.svs.txt"
        proj_df.drop(columns="cbio_project").set_index("Sample_Id").to_csv(
            out, sep="\t", quoting=csv.QUOTE_NONE
        )


if __name__ == "__main__":
    main()
