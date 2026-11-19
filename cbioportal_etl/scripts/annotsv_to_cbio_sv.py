#!/usr/bin/env python3
"""Converts our current DNA annotSV inputs into the cBio SV format.

Focuses on the "split" rows of the Annotation_mode to get per-gene SV data.
Files are processed one at a time and appended to the per-project output.
"""

import argparse
import logging
from collections.abc import Iterator
from pathlib import Path

import pandas as pd
from sv_common import (  # sv_setup_outdir_metadata re-exported for existing importers
    SVWriter,
    configure_logging,
    sv_setup_outdir_metadata,
)

COLUMNS: list[str] = [
    "Sample_Id", "SV_Status", "Site1_Hugo_Symbol", "Site1_Chromosome", "Site1_Contig",
    "Site1_Position", "Site2_Hugo_Symbol", "Site2_Chromosome", "Site2_Contig", "Site2_Position",
    "Site2_Effect_On_Frame", "NCBI_Build", "Class", "Tumor_Split_Read_Count",
    "Tumor_Paired_End_Read_Count", "Event_Info", "Breakpoint_Type", "Annotation", "DNA_Support",
    "SV_Length", "Normal_Paired_End_Read_Count", "Normal_Split_Read_Count", "Comments",
]
RENAME_DICT: dict[str, str] = {
    "AnnotSV_ID": "Event_Info",
    "SV_chrom": "Site1_Chromosome",
    "Tx_start": "Site1_Position",
    "Tx_end": "Site2_Position",
    "Overlapped_tx_length": "SV_Length",
    "CytoBand": "Site1_Contig",
    "Gene_name": "Site1_Hugo_Symbol",
}
SV_TYPE_CLASS: dict[str, str] = {
    "DEL": "Deletion",
    "DUP": "Duplication",
    "INV": "Inversion",
    "INS": "Insertion",
    "BND": "Breakend",
}
# ClinGen annotSV column -> label used in the Annotation column. Add or change entries here;
# the input columns, numeric parsing and annotation text all follow this dict.
CLINGEN_COLUMNS: dict[str, str] = {
    "HI": "Haploinsufficiency",
    "TS": "Triplosensitivity",
}
# every annotSV column used; sample-specific BS ID columns are added per file
STATIC_INPUT_COLS: frozenset[str] = frozenset(
    {"Annotation_mode", "SV_type", "HI", "TS", "ACMG_class", "INFO", "Frameshift", *RENAME_DICT}
)


def make_annotation(df: pd.DataFrame) -> pd.Series:
    """Vectorized ClinGen and ACMG annotation, e.g. 'ClinGen Haploinsufficiency=3;ACMG_class=4'.
 
    Labels and values come from CLINGEN_COLUMNS; a column is skipped for rows where it is blank.
    """
    labels = pd.Series("", index=df.index, dtype=object)
    values = pd.Series("", index=df.index, dtype=object)
    for column, label in CLINGEN_COLUMNS.items():
        col = df[column].astype(object)
        has = col.notna() & (col.astype(str) != "")
        sep = (labels != "").map({True: ",", False: ""})  # comma only after a previous hit
        labels += (sep + label).where(has, "")
        values += (sep + col.where(has, "").astype(str)).where(has, "")
    clingen = ("ClinGen " + labels + "=" + values).where(labels != "", "")
    acmg = df["ACMG_class"].astype(object).fillna("").astype(str).str.removeprefix("full=")
    acmg_part = ("ACMG_class=" + acmg).where(~acmg.isin(["", "NA", "nan"]), "")
    return (clingen + ";" + acmg_part).str.strip(";").astype(object)


def parse_read_support(sample_col: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Parse paired-end and split-read ALT support counts from sample fields (FORMAT/PR:SR).

    Follows form: paired_read_ref,paired_read_alt:split_read_ref,split_read_alt
    """
    if sample_col.str.count(":").gt(1).any():
        too_many_fields = "Unexpected extra ':' fields in sample column"
        raise ValueError(too_many_fields)
    parts = sample_col.str.split(":", n=1, expand=True)
    paired = parts[0].str.split(",").str[1].astype("Int64")
    if parts.shape[1] > 1:
        split = pd.to_numeric(parts[1].str.split(",").str[1], errors="coerce").astype("Int64")
    else:
        split = pd.Series(pd.NA, index=sample_col.index, dtype="Int64")
    return paired, split


def process_annotsv_file(path: Path, metadata: object) -> pd.DataFrame:
    """Read one annotSV file and return it formatted for cBio (COLUMNS only)."""
    affected, reference = str(metadata.affected_bs_id), str(metadata.reference_bs_id)
    wanted = STATIC_INPUT_COLS | {affected, reference}
    try:
        df = pd.read_csv(
            path, sep="\t", keep_default_na=False, na_values=[""], usecols=lambda c: c in wanted
        )
    except ValueError as e:
        err = f"Could not parse {path}: {e}"
        raise RuntimeError(err) from e
    missing = wanted - set(df.columns)
    if missing:
        err = f"Required columns {','.join(sorted(missing))} are missing from {path}"
        raise KeyError(err)

    # keep per-gene ("split") rows only; breakends are currently too complicated to parse
    df = df.loc[(df["Annotation_mode"] == "split") & (df["SV_type"] != "BND")].copy()
    if df.empty:
        return pd.DataFrame(columns=COLUMNS)
    try:
        df["Tumor_Paired_End_Read_Count"], df["Tumor_Split_Read_Count"] = parse_read_support(df[affected])
        df["Normal_Paired_End_Read_Count"], df["Normal_Split_Read_Count"] = parse_read_support(df[reference])
        df = df.drop(columns=[affected, reference])
        df["HI"] = pd.to_numeric(df["HI"], errors="coerce").astype("Int64")
        df["TS"] = pd.to_numeric(df["TS"], errors="coerce").astype("Int64")
        df["Annotation"] = make_annotation(df)
        df["Class"] = df["SV_type"].replace(SV_TYPE_CLASS)
        df["Breakpoint_Type"] = df["INFO"].fillna("").str.contains("IMPRECISE").map(
            {True: "IMPRECISE", False: "PRECISE"}
        )
        df["Site2_Effect_On_Frame"] = df["Frameshift"].fillna("").str.contains("yes").map(
            {True: "Frameshift", False: ""}
        )
        df = df.rename(columns=RENAME_DICT)
        df["Site2_Chromosome"] = df["Site1_Chromosome"]
        df["Site2_Hugo_Symbol"] = df["Site1_Hugo_Symbol"]
        df["Site2_Contig"] = df["Site1_Contig"]
        # split contigs that have a delimiter signifying two values
        mask = df["Site1_Contig"].str.contains("-", na=False)
        if mask.any():
            split_vals = df.loc[mask, "Site1_Contig"].str.split("-", n=1, expand=True)
            df.loc[mask, "Site1_Contig"] = split_vals[0]
            df.loc[mask, "Site2_Contig"] = split_vals[1]
        df = df.astype({"Site1_Position": int, "Site2_Position": int, "SV_Length": "Int64"})
    except (KeyError, ValueError, TypeError, AttributeError) as e:
        err = f"Error parsing read support / annotation info in {path}"
        raise RuntimeError(err) from e
    return df.assign(
        Sample_Id=metadata.cbio_sample_name,
        SV_Status="SOMATIC",
        NCBI_Build="GRCh38",
        DNA_Support="Yes",
        Comments="Manta",
    )[COLUMNS]


def iter_annotsv_chunks(sv_results: str, sv_metadata: pd.DataFrame) -> Iterator[tuple[str, pd.DataFrame]]:
    """Yield (cbio_project, formatted_df) for each annotSV file in the metadata."""
    logger = logging.getLogger(__name__)
    logger.info("Processing %d annotSV files from %s", len(sv_metadata), sv_results)
    for i, metadata in enumerate(sv_metadata.itertuples(), start=1):
        yield metadata.cbio_project, process_annotsv_file(Path(sv_results) / metadata.file_name, metadata)
        if i % 25 == 0:
            logger.info("Processed %d of %d annotSV files", i, len(sv_metadata))


def main() -> None:
    """Convert annotSV results to cBio SV format."""
    parser = argparse.ArgumentParser(description="Convert annotSV results to cbio format.")
    parser.add_argument("-t", "--table", required=True,
                        help="Table with cbio project, kf bs ids, cbio IDs, and file names")
    parser.add_argument("-l", "--link-input-dir", default="annotSV_results/",
                        help="DNA SV directory name to store annotSV symlinks")
    parser.add_argument("-o", "--out-dir", default="merged_sv/",
                        help="Result output dir. Default is merged_sv")
    parser.add_argument("-s", "--skip-symlink", action="store_true",
                        help="Flag to skip sym links if repeating the process with the same input files")
    args = parser.parse_args()
    configure_logging()

    meta = sv_setup_outdir_metadata(args.link_input_dir, args.out_dir, args.table, "sv", args.skip_symlink)
    writer = SVWriter(args.out_dir, COLUMNS)
    for project, df in iter_annotsv_chunks(args.link_input_dir, meta):
        writer.write(project, df)
    writer.finish(meta.cbio_project.unique())


if __name__ == "__main__":
    main()
