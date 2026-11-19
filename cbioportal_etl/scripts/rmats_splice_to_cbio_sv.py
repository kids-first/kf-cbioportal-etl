#!/usr/bin/env python3
"""Converts rMATs inputs into the cBio SV format.

Multiple files per sample comprise different splice types.
Files are processed one at a time and appended to the per-project output.
"""

import argparse
import logging
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd
from sv_common import SVWriter, configure_logging, sv_setup_outdir_metadata

DEFAULT_COMMENT = "rMATS v4.1.2, GENCODE v39"
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
RENAME_DICT = {
    "strand": "RMATS.strand",
    "IJC_SAMPLE_1": "RMATS.IJC_SAMPLE_1",
    "SJC_SAMPLE_1": "RMATS.SJC_SAMPLE_1",
    "IncLevel1": "RMATS.IncLevel1",
}
COLUMNS: list[str] = [
    "Sample_Id", "SV_Status", "Site1_Hugo_Symbol", "Site1_Chromosome", "Site1_Position",
    "Site1_Region", "Site2_Hugo_Symbol", "Site2_Chromosome", "Site2_Position", "Site2_Region",
    "Site2_Effect_On_Frame", "NCBI_Build", "Class", "Tumor_Read_Count", "Event_Info", "Annotation",
    "RNA_Support", "SV_Length", "Comments", "RMATS.strand", "RMATS.IJC_SAMPLE_1",
    "RMATS.SJC_SAMPLE_1", "RMATS.IncLevel1",
]


def make_annotation_event_info(
    df: pd.DataFrame, event_type: str, start_col: str, end_col: str
) -> tuple[pd.Series, pd.Series]:
    """Build Annotation and Event_Info series from the raw rMATS columns."""
    annotation = pd.Series("", index=df.index, dtype=object)
    for col in ANNOTATION_COLS[event_type]:
        annotation += f"{col}: " + df[col].astype(str) + ";"
    event_info = (
        df["geneSymbol"].astype(str) + f"_{event_type}_" + df["chr"].astype(str) + ":"
        + df[start_col].astype(str) + "-" + df[end_col].astype(str)
    )
    return annotation, event_info


def process_rmats_file(path: Path, metadata: object, event_type: str, comment: str) -> pd.DataFrame:
    """Read one rMATS file and return it formatted for cBio (COLUMNS only)."""
    start_col, end_col = POSITION_COLS[event_type]
    usecols = [*BASE_COLS, start_col, end_col, *ANNOTATION_COLS[event_type]]
    try:
        df = pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""], usecols=usecols)
    except ValueError as e:
        err = f"Could not parse {path}: {e}"
        raise RuntimeError(err) from e

    start, end = df[start_col], df[end_col]
    sv_length = end - start
    annotation, event_info = make_annotation_event_info(df, event_type, start_col, end_col)
    region = "INTRON" if event_type == "RETAINED_INTRON" else "EXON"
    if event_type == "RETAINED_INTRON":
        frame = "NA"
    else:
        frame = np.where(sv_length % 3 == 0, "In_frame", "Frameshift")
    return df.rename(columns=RENAME_DICT).assign(
        Sample_Id=metadata.cbio_sample_name,
        SV_Status="SOMATIC",
        Site1_Hugo_Symbol=df["geneSymbol"],
        Site1_Chromosome=df["chr"],
        Site1_Position=start,
        Site1_Region=region,
        Site2_Hugo_Symbol=df["geneSymbol"],
        Site2_Chromosome=df["chr"],
        Site2_Position=end,
        Site2_Region=region,
        Site2_Effect_On_Frame=frame,
        NCBI_Build="GRCh38",
        Class=event_type,
        Tumor_Read_Count=df["IJC_SAMPLE_1"] + df["SJC_SAMPLE_1"],
        Event_Info=event_info,
        Annotation=annotation,
        RNA_Support="Yes",
        SV_Length=sv_length,
        Comments=comment,
    )[COLUMNS]


def iter_rmats_chunks(
    splice_results: str, splice_metadata: pd.DataFrame, comment: str = DEFAULT_COMMENT
) -> Iterator[tuple[str, pd.DataFrame]]:
    """Yield (cbio_project, formatted_df) for each rMATS file in the metadata."""
    logger = logging.getLogger(__name__)
    logger.info("Processing %d rmats files from %s", len(splice_metadata), splice_results)
    for i, metadata in enumerate(splice_metadata.itertuples(), start=1):
        event_type = next((code for token, code in EVENT_TYPES.items() if token in metadata.file_name), None)
        if event_type is None:
            err = f"Could not find valid event type for {metadata.file_name}"
            raise ValueError(err)
        path = Path(splice_results) / metadata.file_name
        yield metadata.cbio_project, process_rmats_file(path, metadata, event_type, comment)
        if i % 25 == 0:
            logger.info("Processed %d of %d rmats files", i, len(splice_metadata))


def main() -> None:
    """Convert rmats splice results to cBio SV format."""
    parser = argparse.ArgumentParser(description="Convert rmats splice results to cbio format.")
    parser.add_argument("-t", "--table", required=True,
                        help="Table with cbio project, kf bs ids, cbio IDs, and file names")
    parser.add_argument("-l", "--link-input-dir", default="rmats_results/",
                        help="rmats splice directory name to store symlinks")
    parser.add_argument("-o", "--out-dir", default="merged_sv/",
                        help="Result output dir. Default is merged_sv")
    parser.add_argument("-s", "--skip-symlink", action="store_true",
                        help="Flag to skip sym links if repeating the process with the same input files")
    parser.add_argument("--comment", default=DEFAULT_COMMENT,
                        help="rMATS and GENCODE version information for the Comments column")
    args = parser.parse_args()
    configure_logging()

    meta = sv_setup_outdir_metadata(args.link_input_dir, args.out_dir, args.table, "splice", args.skip_symlink)
    writer = SVWriter(args.out_dir, COLUMNS)
    for project, df in iter_rmats_chunks(args.link_input_dir, meta, args.comment):
        writer.write(project, df)
    writer.finish(meta.cbio_project.unique())


if __name__ == "__main__":
    main()
