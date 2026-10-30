#!/usr/bin/env python3
"""Converts our current DNA annotSV inputs into the cBio SV format.

Focuses on the "split" rows of the Annotation_mode to get per-gene SV data
"""

import argparse
import csv
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def sv_setup_outdir_metadata(link_input_dir: str, out_dir: str, table: str, skip: bool) -> pd.DataFrame:
    """Create output dir and link input dir if they don't exist, then subset metadata for DNA SV files.

    Output dir will store the merged SV files, and link_input_dir will store symlinks to the annotSV results for processing.

    Args:
        link_input_dir: Directory to store symlinks to annotSV results
        out_dir: Directory to store merged SV results
        table: Table with cbio project, kf bs ids, cbio IDs, and file names
        skip: Flag to skip symlinking if repeating the process with the same input files
    Returns:
        Pandas dataframe with DNA SV-specific metadata

    """
    logger = logging.getLogger(__name__)
    # create output dir and link input dir if they don't exist
    Path(link_input_dir).mkdir(exist_ok=True, parents=True)
    Path(out_dir).mkdir(exist_ok=True, parents=True)

    # deal only with SV metadata
    ext = "sv"
    # ensure sample name is imported as str
    all_file_meta: pd.DataFrame = pd.read_csv(table, sep="\t", dtype={"cbio_sample_name": str})
    dna_sv_subset: pd.DataFrame = all_file_meta.loc[all_file_meta["etl_file_type"] == ext]
    # Create symlinks to mafs in one place for ease of processing
    sv_dirs: list[str] = dna_sv_subset.file_type.unique().tolist()
    logger.info(f"{'Skipping ' if skip else ''}Symlinking SV files from {','.join(sv_dirs)} to {link_input_dir}")
    if skip:
        sv_dirs = []
    sym_errs: int = 0
    for dirname in sv_dirs:
        if not Path(dirname).exists():
            input_err= f"Input SV dir {dirname} does not exist. Check if files were downloaded to the correct location"
            raise ValueError(input_err)
        abs_path: Path = Path(dirname).resolve()
        try:
            for fname in Path(dirname).iterdir():
                src: str = str(abs_path / fname)
                dest: str = str(Path(link_input_dir) / fname)
                Path(dest).symlink_to(src)
        except Exception as e:
            logger.exception("Could not sym link %s in %s", fname, dirname)
            sym_errs += 1
    # If symlink errors, stop here as data will be incomplete
    if sym_errs:
        sym_err_msg = "Could not sym link %d files, exiting!", sym_errs
        raise FileNotFoundError(sym_err_msg)
    return dna_sv_subset


def rename_and_format_cols(concat_sv_df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns and parse data."""
    logger = logging.getLogger(__name__)
    logger.info("Renaming and formatting annotSV columns for cBio SV output")
    desired: list[str] = [
        "Sample_Id",
        "SV_Status",
        "Site1_Hugo_Symbol",
        "Site1_Chromosome",
        "Site1_Contig",
        "Site1_Position",
        "Site2_Hugo_Symbol",
        "Site2_Chromosome",
        "Site2_Contig",
        "Site2_Position",
        "Site2_Effect_On_Frame",
        "NCBI_Build",
        "Class",
        "Tumor_Split_Read_Count",
        "Tumor_Paired_End_Read_Count",
        "Event_Info",
        "Breakpoint_Type",
        "Annotation",
        "DNA_Support",
        "SV_Length",
        "Normal_Paired_End_Read_Count",
        "Normal_Split_Read_Count",
        "Comments",
    ]
    rename_dict: dict[str, str] = {
        "AnnotSV_ID": "Event_Info",
        "SV_chrom": "Site1_Chromosome",
        "Tx_start": "Site1_Position",
        "Tx_end": "Site2_Position",
        "Overlapped_tx_length": "SV_Length",
        "CytoBand": "Site1_Contig",
        "Gene_name": "Site1_Hugo_Symbol",
    }
    sv_type_class_dict: dict[str, str] = {
        "DEL": "Deletion",
        "DUP": "Duplication",
        "INV": "Inversion",
        "INS": "Insertion",
        "BND": "Breakend",
    }
    concat_sv_df["Class"] = concat_sv_df["SV_type"].replace(sv_type_class_dict)
    concat_sv_df["Breakpoint_Type"] = np.where(
        concat_sv_df["INFO"].fillna("").str.contains("IMPRECISE"),
        "IMPRECISE",
        "PRECISE",
    )
    concat_sv_df["Site2_Effect_On_Frame"] = np.where(
        concat_sv_df["Frameshift"].fillna("").str.contains("yes"),
        "Frameshift",
        "",
    )
    concat_sv_df.rename(columns=rename_dict, inplace=True)
    concat_sv_df["Site2_Chromosome"] = concat_sv_df["Site1_Chromosome"]
    concat_sv_df["Site2_Hugo_Symbol"] = concat_sv_df["Site1_Hugo_Symbol"]
    concat_sv_df["Site2_Contig"] = concat_sv_df["Site1_Contig"]
    # split contigs that have a delimeter signifying two values
    mask = concat_sv_df["Site1_Contig"].str.contains("-", na=False)
    split_vals = concat_sv_df.loc[mask, "Site1_Contig"].str.split("-", n=1, expand=True)
    concat_sv_df.loc[mask, "Site1_Contig"] = split_vals[0]
    concat_sv_df.loc[mask, "Site2_Contig"] = split_vals[1]

    concat_sv_df = concat_sv_df.astype({"Site1_Position": int, "Site2_Position": int, "SV_Length": "Int64"})
    # set some defaults
    concat_sv_df["SV_Status"] = "SOMATIC"
    concat_sv_df["NCBI_Build"] = "GRCh38"
    concat_sv_df["DNA_Support"] = "Yes"
    concat_sv_df["Comments"] = "Manta"
    return concat_sv_df[desired]


def make_annotation(row: pd.Series) -> str:
    """Parse ClinGen HI, TS and ACMG_Class to make a single annotation."""
    parts = []

    clingen_columns = {
        "HI": "Haploinsufficiency",
        "TS": "Triplosensitivity",
    }
    clingen_hits = [
        (label, row[column])
        for column, label in clingen_columns.items()
        if pd.notna(row[column]) and row[column] != ""
    ]
    if clingen_hits:
        labels, values = zip(*clingen_hits, strict=True)
        parts.append(f"ClinGen {','.join(labels)}={','.join(map(str, values))}")

    acmg = str(row["ACMG_class"]).removeprefix("full=")
    if acmg and acmg not in {"NA", "nan"}:
        parts.append(f"ACMG_class={acmg}")

    return ";".join(parts)


def init_cbio_sv_df(sv_results: str, sv_metadata: pd.DataFrame) -> pd.DataFrame:
    """Use data frame subset on DNA SV files to find and merge result files.

    Args:
        sv_results: annotSV results dir
        sv_metadata: Pandas dataframe with DNA SV-specific metadata
    Returns:
        Collapsed and formatted dataframe, list of desired fields

    """
    logger = logging.getLogger(__name__)
    logger.info("Processing %d annotSV files from %s", len(sv_metadata), sv_results)
    frame_list = []
    for metadata in sv_metadata.itertuples():
        # concat annot SV file, rename Sample Column according to cBio name
        try:
            fname = f"{sv_results}/{metadata.file_name}"
            ann_file: pd.DataFrame = pd.read_csv(
                fname,
                sep="\t", keep_default_na=False, na_values=[""]
            )
        except Exception as e:
            file_process_err = "Could not open files %s in %s", metadata.file_name, sv_results
            raise FileNotFoundError(file_process_err) from e
        # validate keys/columns are there before proceeding
        required_cols: set[str] = {"Annotation_mode", str(metadata.affected_bs_id),
                                   str(metadata.reference_bs_id), "SV_type", "HI", "TS" }
        if not required_cols.issubset(set(ann_file.columns)):
            key_err = f"Required columns {",".join(set(required_cols).difference(ann_file.columns))} are missing dataframe"
            raise KeyError(key_err)
        # drop entries where Annotation_mode is not "split" (i.e. per-gene)
        # and not breakend (currently too complicated to parse)
        ann_file: pd.DataFrame = ann_file.loc[(ann_file["Annotation_mode"] == "split") & (ann_file["SV_type"] != "BND")]
        ann_file = ann_file.assign(Sample_Id=metadata.cbio_sample_name)
        try:
            # at this step, get the tumor BS ID
            # parse FORMAT for PR and SR and assign values cBio columns
            # drop the tumor and normal BS ID cols
            affected_read_support = ann_file[metadata.affected_bs_id].str.split(":", expand=True)
            reference_read_support = ann_file[metadata.reference_bs_id].str.split(":", expand=True)
            ann_file["Tumor_Paired_End_Read_Count"] = affected_read_support[0].str.split(",").str[1].astype("Int64")
            ann_file["Normal_Paired_End_Read_Count"] = reference_read_support[0].str.split(",").str[1].astype("Int64")
            ann_file["Tumor_Split_Read_Count"] = pd.Series(pd.NA, index=ann_file.index, dtype="Int64")
            ann_file["Normal_Split_Read_Count"] = pd.Series(pd.NA, index=ann_file.index, dtype="Int64")

            if len(affected_read_support.columns) > 1:
                ann_file["Tumor_Split_Read_Count"] = pd.to_numeric(
                    affected_read_support[1].str.split(",").str[1],
                    errors="coerce",
                ).astype("Int64")

                ann_file["Normal_Split_Read_Count"] = pd.to_numeric(
                    reference_read_support[1].str.split(",").str[1],
                    errors="coerce",
                ).astype("Int64")
            ann_file = ann_file.drop(columns=[metadata.affected_bs_id, metadata.reference_bs_id])
        except Exception as e:
            ct_parse_err = f"Error parsing sample FORMAT columns for read support in {fname}"
            raise RuntimeError(ct_parse_err) from e
        try:
            # Create Annotation column parsing HI, TS, ACMG_class
            ann_file["HI"] = pd.to_numeric(ann_file["HI"], errors="coerce").astype("Int64")
            ann_file["TS"] = pd.to_numeric(ann_file["TS"], errors="coerce").astype("Int64")
            ann_file["Annotation"] = ann_file.apply(make_annotation, axis=1)
            frame_list.append(ann_file)
        except Exception as e:
            annot_err = f"Error parsing ClinGen and ACMG annotation info in {fname}"
            raise RuntimeError(annot_err) from e
    concat_frame: pd.DataFrame = pd.concat(frame_list)
    del frame_list
    return concat_frame


def main() -> None:
    """Convert annotSV results to cBio SV format."""
    parser = argparse.ArgumentParser(
        description="Convert annotSV results to cbio format.",
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
    parser.add_argument(
        "-s",
        "--skip-symlink",
        action="store_true",
        help="Flag to skip sym links if repeating the process with the same input files",
        required=False,
    )
    args = parser.parse_args()

    dna_sv_subset: pd.DataFrame = sv_setup_outdir_metadata(args.link_input_dir, args.out_dir, args.table, args.skip_symlink)
    project_list: np.ndarray = dna_sv_subset.cbio_project.unique()
    cbio_sv_df: pd.DataFrame = init_cbio_sv_df(args.link_input_dir, dna_sv_subset)
    cbio_sv_df = rename_and_format_cols(concat_sv_df=cbio_sv_df)


    for project in project_list:
        cbio_sv_fname: Path = Path(args.out_dir) / (project + ".svs.txt")
        cbio_sv_df.set_index("Sample_Id", inplace=True)
        cbio_sv_df.to_csv(cbio_sv_fname, sep="\t", mode="w", index=True, quoting=csv.QUOTE_NONE)


if __name__ == "__main__":
    main()
