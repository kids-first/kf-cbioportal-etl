#!/usr/bin/env python3
"""Converts our current fusion inputs into the cBio SV format.

For repeat rows (same breakpoint in a sample, different callers),
ARRIBA annot is used with ceiling of mean for counts used.
Samples (annoFuse) or batches of samples (openX/dgd) are processed one at a time.
"""

import argparse
import logging
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd
from sv_common import SVWriter, configure_logging

SAMPLES_PER_BATCH = 100  # openX/dgd: samples processed per chunk
KEY_COLS: list[str] = [
    "Sample", "Gene1A", "LeftBreakpoint", "Gene1B", "RightBreakpoint", "FusionName", "Fusion_anno",
]
SORT_COLS: list[str] = ["Sample", "LeftBreakpoint", "RightBreakpoint", "Caller"]
DESIRED_INPUT: list[str] = [
    "Sample", "Gene1A", "LeftBreakpoint", "Gene1B", "RightBreakpoint", "Fusion_Type",
    "JunctionReadCount", "SpanningFragCount", "annots", "FusionName", "Fusion_anno", "Caller",
]
# every input column that is used anywhere; everything else is never read from disk
INPUT_COLS = frozenset({*DESIRED_INPUT, "CalledBy", "Gene1A_anno", "Gene1B_anno", "Confidence"})
RENAME_DICT: dict[str, str] = {
    "Sample": "Sample_Id",
    "Gene1A": "Site1_Hugo_Symbol",
    "Gene1B": "Site2_Hugo_Symbol",
    "Fusion_Type": "Site2_Effect_On_Frame",
    "JunctionReadCount": "Tumor_Read_Count",
    "SpanningFragCount": "Tumor_Split_Read_Count",
    "annots": "Annotation",
    "FusionName": "Event_Info",
    "Fusion_anno": "External_Annotation",
    "Caller": "Comments",
}
COLUMNS: list[str] = [
    "Sample_Id", "SV_Status", "Site1_Hugo_Symbol", "Site1_Chromosome", "Site1_Position",
    "Site2_Hugo_Symbol", "Site2_Chromosome", "Site2_Position", "Site2_Effect_On_Frame",
    "NCBI_Build", "Tumor_Read_Count", "Tumor_Split_Read_Count", "Annotation", "RNA_Support",
    "Connection_Type", "Event_Info", "Class", "External_Annotation", "Comments",
]


def fusion_setup_outdir_metadata(mode: str, out_dir: str, table: str) -> pd.DataFrame:
    """Subset metadata based on run mode (openX/dgd/kfprod) and create the output dir."""
    Path(out_dir).mkdir(exist_ok=True, parents=True)
    r_ext = {"openX": "rsem", "dgd": "DGD_FUSION"}.get(mode, "fusion")
    # ensure sample name is imported as str
    all_file_meta = pd.read_csv(table, sep="\t", dtype={"cbio_sample_name": str})
    rna_subset = all_file_meta.loc[all_file_meta["etl_file_type"] == r_ext]
    if rna_subset.empty and r_ext != "rsem":
        r_ext = "rsem"
        rna_subset = all_file_meta.loc[all_file_meta["etl_file_type"] == r_ext]
    if rna_subset.empty:
        err = f"No metadata obtained using category {r_ext}"
        raise KeyError(err)
    return rna_subset.reset_index(drop=True)


def _collapse_groups(
    df: pd.DataFrame, join_cols: dict[str, str], ceil_mean_cols: tuple[str, ...] = ()
) -> pd.DataFrame:
    """Collapse rows sharing KEY_COLS into their first row (after sorting by caller).

    Args:
        df: Fusion rows
        join_cols: column -> separator; the unique values in each group are joined
        ceil_mean_cols: columns replaced by the ceiling of the group mean
    """
    df = df.sort_values(SORT_COLS, ignore_index=True)
    gid = df.groupby(KEY_COLS, sort=False, dropna=False).ngroup()
    first = ~gid.duplicated()
    out = df.loc[first].copy()
    first_gid = gid[first].to_numpy()
    grouped = df.groupby(gid)
    for col, sep in join_cols.items():
        joined = grouped[col].agg(lambda s, sep=sep: sep.join(dict.fromkeys(s.dropna().astype(str))))
        out[col] = joined.reindex(first_gid).to_numpy()
    for col in ceil_mean_cols:
        out[col] = np.ceil(grouped[col].mean().reindex(first_gid).to_numpy()).astype(int)
    return out


def collapse_and_format(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse the same call from different callers; counts become ceiling of the mean."""
    counts = ["JunctionReadCount", "SpanningFragCount"]
    out = _collapse_groups(df, {"Caller": ","}, tuple(counts))
    out["Caller"] = out["Caller"].str.upper()
    carry = [c for c in ("annots", "Fusion_Type") if c in out.columns]
    return out[[*KEY_COLS, *counts, *carry, "Caller"]]


def openx_annot_bug(df: pd.DataFrame) -> pd.DataFrame:
    """Merge same fusion calls with annotation split across multiple lines."""
    out = _collapse_groups(df, {"Caller": ",", "Gene1A_anno": ", ", "Gene1B_anno": ", "})
    out["Caller"] = out["Caller"].str.upper()
    return out


def _unique_join(s: pd.Series) -> pd.Series:
    return s.fillna("").map(lambda v: ",".join(dict.fromkeys(v.split(", "))))


def build_annots(df: pd.DataFrame) -> pd.Series:
    """Prefix annots with the (de-duplicated) Gene1A/Gene1B annotations."""
    return (
        "Gene1: " + _unique_join(df["Gene1A_anno"])
        + "; Gene2: " + _unique_join(df["Gene1B_anno"])
        + ";" + df["annots"]
    )


def final_order_cols(df: pd.DataFrame, mode: str) -> pd.DataFrame:
    """Rename columns, set defaults, and order for cBio."""
    df = df.rename(columns=RENAME_DICT).assign(
        Class="FUSION", SV_Status="SOMATIC", NCBI_Build="GRCh38", Connection_Type="5to3",
        RNA_Support="Yes", Site1_Chromosome="", Site1_Position="", Site2_Chromosome="", Site2_Position="",
    )
    if mode != "dgd":
        df[["Site1_Chromosome", "Site1_Position"]] = df["LeftBreakpoint"].str.split(":", expand=True)
        df[["Site2_Chromosome", "Site2_Position"]] = df["RightBreakpoint"].str.split(":", expand=True)
    # NOTE: preserved from the original script: Fusion_Type is always blanked here. If that is
    # unintended, make this `if mode == "dgd"` and keep upper()/replace("-", "_") for the rest.
    df["Site2_Effect_On_Frame"] = ""
    return df.reindex(columns=COLUMNS)


def _finalize(df: pd.DataFrame, mode: str, *, drop_low: bool) -> pd.DataFrame:
    """Filter, annotate, collapse, and format one chunk of fusion rows."""
    df = df.loc[df["Confidence"] != "low"] if drop_low else openx_annot_bug(df)
    if df.empty:
        return pd.DataFrame(columns=COLUMNS)
    df = df.assign(annots=build_annots(df))
    df = df[[c for c in DESIRED_INPUT if c in df.columns]]
    # only collapse counts if read counts are there
    if drop_low or "JunctionReadCount" in df.columns:
        df = collapse_and_format(df)
    return final_order_cols(df, mode)


def _read(path: str | Path) -> pd.DataFrame:
    try:
        return pd.read_csv(
            path, sep="\t", keep_default_na=False, na_values=[""], usecols=lambda c: c in INPUT_COLS
        )
    except ValueError as e:
        err = f"Could not parse {path}: {e}"
        raise RuntimeError(err) from e


def _iter_openx_chunks(path: str, mode: str, rna_metadata: pd.DataFrame) -> Iterator[tuple[str, pd.DataFrame]]:
    openx = _read(path)
    # Merge so that sample names can be cBio names; rows without metadata are dropped
    meta = rna_metadata[["affected_bs_id", "cbio_sample_name", "cbio_project"]].drop_duplicates()
    merged = openx.merge(meta, left_on="Sample", right_on="affected_bs_id", how="inner")
    del openx
    merged["Sample"] = merged["cbio_sample_name"]
    merged = merged.drop(columns=["affected_bs_id", "cbio_sample_name"])
    # openPBTA...and maybe open pedcan uses this
    if "CalledBy" in merged.columns:
        merged = merged.rename(columns={"CalledBy": "Caller"})
    elif mode == "dgd":
        merged["Caller"] = "Archer"
    for project, proj_df in merged.groupby("cbio_project", sort=False):
        # process whole samples at a time so collapsing stays correct and memory stays bounded
        batch_ids = pd.factorize(proj_df["Sample"])[0] // SAMPLES_PER_BATCH
        for _, batch in proj_df.groupby(batch_ids, sort=False):
            yield project, _finalize(batch.drop(columns="cbio_project"), mode, drop_low=False)


def _iter_annofuse_chunks(
    results_dir: str, mode: str, rna_metadata: pd.DataFrame
) -> Iterator[tuple[str, pd.DataFrame]]:
    logger = logging.getLogger(__name__)
    # all files of one sample are collapsed together, as the original concat-then-collapse did
    groups = rna_metadata.groupby(["cbio_project", "cbio_sample_name"], sort=False)
    for i, ((project, sample), grp) in enumerate(groups, start=1):
        frames = [_read(Path(results_dir) / name) for name in grp["file_name"]]
        df = pd.concat(frames, ignore_index=True).assign(Sample=sample)
        yield project, _finalize(df, mode, drop_low=True)
        if i % 25 == 0:
            logger.info("Processed %d of %d fusion samples", i, groups.ngroups)


def iter_fusion_chunks(
    fusion_results: str, mode: str, rna_metadata: pd.DataFrame
) -> Iterator[tuple[str, pd.DataFrame]]:
    """Yield (cbio_project, formatted_df) chunks from annoFuse files (kfprod) or a merged openX/dgd file."""
    if mode in {"openX", "dgd"}:
        yield from _iter_openx_chunks(fusion_results, mode, rna_metadata)
    else:
        yield from _iter_annofuse_chunks(fusion_results, mode, rna_metadata)


def main() -> None:
    """Convert openX fusion or list of annoFuse files to cBio format."""
    parser = argparse.ArgumentParser(
        description="Convert openX fusion table OR list of annofuse files to cbio format."
    )
    parser.add_argument("-t", "--table", required=True,
                        help="Table with cbio project, kf bs ids, cbio IDs, and file names")
    parser.add_argument("-f", "--fusion-results", required=True,
                        help="annoFuse results dir OR openX merged fusion file")
    parser.add_argument("-o", "--out-dir", default="merged_sv/",
                        help="Result output dir. Default is merged_sv")
    parser.add_argument("-m", "--mode", choices=["openX", "kfprod", "dgd"], required=True,
                        help="describe source")
    parser.add_argument("-a", "--append", action="store_true",
                        help="Append to existing project files, reusing their header")
    args = parser.parse_args()
    configure_logging()

    meta = fusion_setup_outdir_metadata(args.mode, args.out_dir, args.table)
    writer = SVWriter(args.out_dir, COLUMNS, append=args.append)
    for project, df in iter_fusion_chunks(args.fusion_results, args.mode, meta):
        writer.write(project, df)
    writer.finish(meta.cbio_project.unique())


if __name__ == "__main__":
    main()
