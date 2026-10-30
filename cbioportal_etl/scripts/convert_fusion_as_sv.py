#!/usr/bin/env python3
"""Converts our current fusion inputs into the cBio SV format.

For repeat rows (same breakpoint in a sample, different callers),
ARRIBA annot is used with ceiling of mean for counts used
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

def fusion_setup_outdir_metadata(mode: str, out_dir: str, table: str) -> pd.DataFrame:
    """Subset metadata based on run mode.

    Args:
        mode: openX if from am openPBTA/pedcan project, standard piplines, or DGD
        out_dir: Name of output directorry to make
        table: cBio formatted metadata table

    Returns:
        Subsetted dataframe

    """
    Path(out_dir).mkdir(exist_ok=True, parents=True)

    # deal only with RNA metadata
    r_ext: str = "fusion"
    if mode == "openX":
        r_ext = "rsem"
    elif mode == "dgd":
        r_ext = "DGD_FUSION"
    # ensure sample name is imported as str
    all_file_meta: pd.DataFrame = pd.read_csv(table, sep="\t", dtype={"cbio_sample_name": str})
    rna_subset =  all_file_meta.loc[all_file_meta["etl_file_type"] == r_ext]
    if rna_subset is None and r_ext != "rsem":
        r_ext = "rsem"
        rna_subset = all_file_meta.loc[all_file_meta["etl_file_type"] == r_ext]
    else:
        no_meta_err = f"No metadata obtained using category {r_ext}"
        raise KeyError(no_meta_err)
    # reset index so that references work later while iterating
    return rna_subset.reset_index(drop=True)


def final_order_cols(fusion_data_collapsed: pd.DataFrame, present_cols: list[str], mode:str) -> pd.DataFrame:
    """Rename columns set defaults."""
    # Get relevant columns
    fusion_data_collapsed.set_index("Sample", inplace=True)
    fusion_data_collapsed.reset_index(inplace=True)
    # openX and annoFuse have different, cols
    rename_dict: dict[str, str] = {
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
    fusion_data_collapsed.rename(columns=rename_dict, inplace=True)
    # Fill in some defaults
    fusion_data_collapsed["Class"] = "FUSION"
    fusion_data_collapsed["SV_Status"] = "SOMATIC"
    fusion_data_collapsed["NCBI_Build"] = "GRCh38"
    fusion_data_collapsed["Connection_Type"] = "5to3"
    fusion_data_collapsed["RNA_Support"] = "Yes"
    fusion_data_collapsed["Site1_Chromosome"] = ""
    fusion_data_collapsed["Site1_Position"] = ""
    fusion_data_collapsed["Site2_Chromosome"] = ""
    fusion_data_collapsed["Site2_Position"] = ""
    # Split some columns that have 2 cols worth of info
    if mode != "dgd":
        fusion_data_collapsed[["Site1_Chromosome", "Site1_Position"]] = fusion_data_collapsed.LeftBreakpoint.str.split(
            ":", expand=True
        )
        fusion_data_collapsed[["Site2_Chromosome", "Site2_Position"]] = fusion_data_collapsed.RightBreakpoint.str.split(
            ":", expand=True
        )

    # Reformat values to fit needs to be ALL CAPS, replace - with _, remove weird chars
    fusion_data_collapsed["Site2_Effect_On_Frame"] = ""
    if mode != "dgd":
        fusion_data_collapsed["Site2_Effect_On_Frame"] = fusion_data_collapsed["Site2_Effect_On_Frame"].str.upper()
        fusion_data_collapsed["Site2_Effect_On_Frame"] = fusion_data_collapsed["Site2_Effect_On_Frame"].str.replace(
            "-", "_"
        )
    # Reorder table
    order_list: list[str] = [
        "Sample_Id",
        "SV_Status",
        "Site1_Hugo_Symbol",
        "Site1_Entrez_Gene_Id",
        "Site1_Ensembl_Transcript_Id",
        "Site1_Exon",
        "Site1_Chromosome",
        "Site1_Position",
        "Site2_Hugo_Symbol",
        "Site2_Entrez_Gene_Id",
        "Site2_Ensembl_Transcript_Id",
        "Site2_Exon",
        "Site2_Chromosome",
        "Site2_Position",
        "Site2_Effect_On_Frame",
        "NCBI_Build",
        "Tumor_Read_Count",
        "Tumor_Split_Read_Count",
        "Annotation",
        "RNA_Support",
        "Connection_Type",
        "Event_Info",
        "Class",
        "External_Annotation",
        "Comments",
    ]
    # again, ensure only present columns are ordered
    present_order: list[str] = [col for col in order_list if col in fusion_data_collapsed]
    return fusion_data_collapsed[present_order]


def collapse_and_format(fusion_data: pd.DataFrame) -> pd.DataFrame:
    """Collapse same call be different caller into same line.

    Use ceiling of mean values of read and frag counts, then reduce to desired columns
    Args:
        fusion_data: Fusion data dataframe
    Returns:
        Collapsed dataframe
    """
    logger = logging.getLogger(__name__)
    # Sort as a cheat to easily select preferred annot later
    logger.info("Sorting for collapse by caller and counts")
    fusion_data = fusion_data.sort_values(["Sample", "LeftBreakpoint", "RightBreakpoint", "Caller"])
    # Merge rows that have the exact same fusion from different callers - thanks Natasha!
    key_cols: list[str] = [
        "Sample",
        "Gene1A",
        "LeftBreakpoint",
        "Gene1B",
        "RightBreakpoint",
        "FusionName",
        "Fusion_anno",
    ]
    logger.info("Grouping for collapse and counts")
    fusion_data["groupby_key"] = fusion_data.apply(
        lambda row: "\t".join([str(row[col]) for col in key_cols]), axis=1
    )
    # "JunctionReadCount","SpanningFragCount","annots"
    collapsed_list: list = []

    for g in fusion_data.groupby(by="groupby_key"):
        _values, df_group = g
        df_group["Caller"] = ",".join(set(df_group["Caller"].tolist()))
        df_group["JunctionReadCount"] = df_group["JunctionReadCount"].mean()
        df_group["SpanningFragCount"] = df_group["SpanningFragCount"].mean()
        # Go with the ceiling of the mean
        df_group["JunctionReadCount"] = df_group["JunctionReadCount"].apply(np.ceil)
        df_group["SpanningFragCount"] = df_group["SpanningFragCount"].apply(np.ceil)
        collapsed_list.append(df_group[key_cols + [ "JunctionReadCount","SpanningFragCount","annots", "Fusion_Type", "Caller"]].head(1))
    del fusion_data
    fusion_data_collapsed = pd.concat(collapsed_list)
    del collapsed_list
    # Should be int
    fusion_data_collapsed["JunctionReadCount"] = fusion_data_collapsed["JunctionReadCount"].astype(
        int
    )
    fusion_data_collapsed["SpanningFragCount"] = fusion_data_collapsed["SpanningFragCount"].astype(
        int
    )
    fusion_data_collapsed["Caller"] = fusion_data_collapsed["Caller"].str.upper()
    return fusion_data_collapsed


def openx_annot_bug(fusion_data: pd.DataFrame) -> pd.DataFrame:
    """Merge same fusion calls with annotation split across multiple lines."""
    logger = logging.getLogger(__name__)
    # Sort as a cheat to easily select preferred annot later
    logger.info("Sorting for collapse bug fix")
    fusion_data = fusion_data.sort_values(["Sample", "LeftBreakpoint", "RightBreakpoint", "Caller"])
    # Merge rows that have the exact same fusion from different callers - thanks Natasha!
    key_cols: list[str] = [
        "Sample",
        "Gene1A",
        "LeftBreakpoint",
        "Gene1B",
        "RightBreakpoint",
        "FusionName",
        "Fusion_anno",
    ]
    remain_cols: list[str] = list(set(fusion_data.columns.to_list()) - set(key_cols))
    logger.info("Grouping for collapse bug fix")
    fusion_data["groupby_key"] = fusion_data.apply(
        lambda row: "\t".join([str(row[col]) for col in key_cols]),
        axis=1,
    )
    # "JunctionReadCount","SpanningFragCount","annots"
    collapsed_list: list = []
    logger.info("Collapsing annotations for bug fix")
    for g in fusion_data.groupby(by="groupby_key"):
        _values, df_group = g
        df_group["Caller"] = ",".join(set(df_group["Caller"].tolist()))
        df_group["Gene1A_anno"] = ", ".join(set(df_group["Gene1A_anno"].tolist()))
        df_group["Gene1B_anno"] = ", ".join(set(df_group["Gene1B_anno"].tolist()))
        collapsed_list.append(df_group[key_cols + remain_cols].head(1))
    del fusion_data
    fusion_data_collapsed = pd.concat(collapsed_list)
    del collapsed_list
    fusion_data_collapsed["Caller"] = fusion_data_collapsed["Caller"].str.upper()
    logger.info("Bug fix completed")
    return fusion_data_collapsed


def filter_and_format_annots(sample_renamed_df: pd.DataFrame, drop_low: bool) -> pd.DataFrame:
    """Apply a filter to remove entries with ARRIBA confidence "low" when set.

    Also formats annots file to include Gene1A_anno and Gene1B_anno
    Args:
        sample_renamed_df: Dataframe with BS IDs replaced by cBio IDs
        drop_low: Boolean flag to drop low arriba hits
    Returns:
        Filtered and formatted dataframe

    """
    logger = logging.getLogger(__name__)
    if drop_low:
        sample_renamed_df = sample_renamed_df[sample_renamed_df.Confidence != "low"]
    else:
        # not drop low is a OpenX feature, apply repeat annotation bug fix
        logger.info("Applying OpenX annot bug fix")
        sample_renamed_df = openx_annot_bug(sample_renamed_df)
    sample_renamed_df["annots"] = (
        sample_renamed_df.apply(
            lambda row: "Gene1: "
            + ",".join(set(row["Gene1A_anno"].split(", ")))
            + "; Gene2: "
            + ",".join(set(row["Gene1B_anno"].split(", "))),
            axis=1,
        )
        + ";"
        + sample_renamed_df["annots"]
    )
    return sample_renamed_df


def init_cbio_master(
    fusion_results: str, mode: str, rna_metadata: pd.DataFrame) -> pd.DataFrame:
    """Use data frame subset on RNA fusion files to find and merge result files.

    Args:
        fusion_results: annoFuse results dir OR openX merged fusion file depending on the mode
        mode: openX or kfprod, or dgd
        rna_metadata: Pandas dataframe with RNA-specific metadata
    Returns:
        Collapsed and formatted dataframe, list of desired fields

    """
    desired: list[str] = [
        "Sample",
        "Gene1A",
        "LeftBreakpoint",
        "Gene1B",
        "RightBreakpoint",
        "Fusion_Type",
        "JunctionReadCount",
        "SpanningFragCount",
        "annots",
        "FusionName",
        "Fusion_anno",
        "Caller",
    ]
    if mode in {"openX", "dgd"}:
        try:
            openx_data: pd.DataFrame = pd.read_csv(
                fusion_results, sep="\t", keep_default_na=False, na_values=[""]
            )
        except Exception as e:
            missing = "Could not open file %s, make sure this file exists!", fusion_results
            raise FileNotFoundError(missing) from e
        # Merge so that sample names can be cBio names - thanks Natasha!
        merged: pd.DataFrame = pd.merge(
            openx_data,
            rna_metadata[["affected_bs_id", "cbio_sample_name"]],
            left_on="Sample",
            right_on="affected_bs_id",
            how="left",
        )
        merged["Sample"] = merged.apply(lambda row: row["cbio_sample_name"], axis=1)
        # openPBTA...and maybe open pedcan uses this
        if "CalledBy" in merged.columns:
            merged.rename(columns={"CalledBy": "Caller"}, inplace=True)
        elif mode == "dgd":
            merged["Caller"] = "Archer"
        # Also merge existing annotations in Gene1A_anno, Gene1B_anno into annots
        merged = filter_and_format_annots(sample_renamed_df=merged, drop_low=False)
        # OpenX data may not have all annoFuse cols
        present: list[str] = [col for col in desired if col in merged.columns]
        openx_data: pd.DataFrame = merged[present]
        # only if read counts there, collapse
        if "JunctionReadCount" in openx_data.columns:
            openx_data = collapse_and_format(openx_data)
        return final_order_cols(openx_data, present, mode)

    frame_list = []
    try:
        for metadata in rna_metadata.itertuples():
            # concat annofuse file, rename Sample Column according to cBio name
            ann_file: pd.DataFrame = pd.read_csv(
                f"{fusion_results}/{metadata.file_name}", sep="\t", keep_default_na=False, na_values=[""]
            )
            ann_file = ann_file.assign(Sample=metadata.cbio_sample_name)
            frame_list.append(ann_file)
    except Exception as e:
        merge_err = "Could not appaned dataframes and rename sample ID"
        raise RuntimeError(merge_err) from e
    concat_frame: pd.DataFrame = pd.concat(frame_list)
    concat_frame = filter_and_format_annots(sample_renamed_df=concat_frame, drop_low=True)
    concat_frame = collapse_and_format(concat_frame)
    del frame_list
    return final_order_cols(concat_frame, desired, mode)


def main() -> None:
    """Convert openX fusion or list of annoFuse files to cBio format."""
    parser = argparse.ArgumentParser(
        description="Convert openX fusion table OR list of annofuse files to cbio format."
    )
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        dest="table",
        help="Table with cbio project, kf bs ids, cbio IDs, and file names",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--fusion-results",
        action="store",
        dest="fusion_results",
        help="annoFuse results dir OR openX merged fusion file",
        required=True,
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
        "-m",
        "--mode",
        action="store",
        dest="mode",
        help="describe source",
        choices=["openX", "kfprod", "dgd"],
        required=True,
    )
    parser.add_argument(
        "-a",
        "--append",
        action="store_true",
        dest="append",
        help="Flag to append, meaning print to STDOUT and skip header",
        required=False,
    )

    args = parser.parse_args()

    # ext used in pbta vs openpedcan varies
    rna_subset: pd.DataFrame = fusion_setup_outdir_metadata(args.mode, args.out_dir, args.table)

    project_list: np.ndarray = rna_subset.cbio_project.unique()
    cbio_master = init_cbio_master(args.fusion_results, args.mode, rna_subset)

    for project in project_list:
        sub_sample_list = list(
            rna_subset.loc[rna_subset["cbio_project"] == project, "cbio_sample_name"]
        )
        fus_fname: Path = Path(args.out_dir) / (project + ".svs.txt")
        fus_tbl: pd.DataFrame = cbio_master[cbio_master.Sample_Id.isin(sub_sample_list)]
        fus_tbl.fillna("NA", inplace=True)
        if not args.append:
            fus_tbl.set_index("Sample_Id", inplace=True)
            fus_tbl.to_csv(fus_fname, sep="\t", mode="w", index=True, quoting=csv.QUOTE_NONE)
        else:
            # append to existing fusion file, use same header
            existing = pd.read_csv(fus_fname, sep="\t", keep_default_na=False, na_values=[""])
            fus_tbl = fus_tbl[existing.columns]
            fus_tbl.set_index("Sample_Id", inplace=True)
            fus_tbl.to_csv(
                fus_fname, sep="\t", mode="a", index=True, quoting=csv.QUOTE_NONE, header=None
            )


if __name__ == "__main__":
    main()
