#!/usr/bin/env python3
"""Converts our current fusion inputs into the cBio SV format.

For repeat rows (same breakpoint in a sample, different callers),
ARRIBA annot is used with ceiling of mean for counts used
"""

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd


def collapse_and_format(fusion_data: pd.DataFrame) -> pd.DataFrame:
    """Collapse same call be different caller into same line.

    Use ceiling of mean values of read and frag counts
    Args:
        fusion_data: Fusion data dataframe
    Returns:
        Collapsed dataframe
    """
    # Sort as a cheat to easily select preferred annot later
    print("Sorting for collapse by caller and counts", file=sys.stderr)
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
    print("Grouping for collapse and counts", file=sys.stderr)
    fusion_data["groupby_key"] = fusion_data.apply(
        lambda row: "\t".join([str(row[col]) for col in key_cols]), axis=1
    )
    # "JunctionReadCount","SpanningFragCount","annots"
    collapsed_list: list = []

    for g in fusion_data.groupby(by="groupby_key"):
        values, df_group = g
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
    # Sort as a cheat to easily select preferred annot later
    print("Sorting for collapse bug fix", file=sys.stderr)
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
    print("Grouping for collapse bug fix", file=sys.stderr)
    fusion_data["groupby_key"] = fusion_data.apply(
        lambda row: "\t".join([str(row[col]) for col in key_cols]),
        axis=1,
    )
    # "JunctionReadCount","SpanningFragCount","annots"
    collapsed_list: list = []
    print("Collapsing annotations for bug fix", file=sys.stderr)
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
    print("Bug fix completed", file=sys.stderr)
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
    if drop_low:
        sample_renamed_df = sample_renamed_df[sample_renamed_df.Confidence != "low"]
    else:
        # not drop low is a OpenX feature, apply repeat annotation bug fix
        print("Applying OpenX annot bug fix", file=sys.stderr)
        sample_renamed_df = openx_annot_bug(sample_renamed_df)
    sample_renamed_df["annots"] = (
        sample_renamed_df.apply(
            lambda row: "Gene1: "
            + ",".join(set(list(row["Gene1A_anno"].split(", "))))
            + "; Gene2: "
            + ",".join(set(list(row["Gene1B_anno"].split(", ")))),
            axis=1,
        )
        + ";"
        + sample_renamed_df["annots"]
    )
    return sample_renamed_df


def init_cbio_master(
    fusion_results: str, mode: str, rna_metadata: pd.DataFrame
) -> tuple[pd.DataFrame, list[str]]:
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
        openx_data: pd.DataFrame = pd.read_csv(
            fusion_results, sep="\t", keep_default_na=False, na_values=[""]
        )
        # Merge so that sample names can be cBio names - thanks Natasha!
        merged: pd.DataFrame = pd.merge(
            openx_data,
            rna_metadata[["T_CL_BS_ID", "Cbio_Tumor_Name"]],
            left_on="Sample",
            right_on="T_CL_BS_ID",
            how="left",
        )
        merged["Sample"] = merged.apply(lambda row: row["Cbio_Tumor_Name"], axis=1)
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
        return openx_data, present
    else:
        flist: pd.Series[str] = rna_metadata.File_Name
        frame_list = []
        try:
            for i in range(0, len(flist), 1):
                # concat annofuse file, rename Sample Column according to cBio name
                ann_file: pd.DataFrame = pd.read_csv(
                    f"{fusion_results}/{flist[i]}", sep="\t", keep_default_na=False, na_values=[""]
                )
                ann_file = ann_file.assign(Sample=rna_metadata.at[i, "Cbio_Tumor_Name"])
                frame_list.append(ann_file)
        except Exception as e:
            print(f"{e}", file=sys.stderr)
            sys.exit(1)
        concat_frame: pd.DataFrame = pd.concat(frame_list)
        concat_frame = filter_and_format_annots(sample_renamed_df=concat_frame, drop_low=True)
        del frame_list
        fusion_data: pd.DataFrame = concat_frame[desired]
        del concat_frame
        return collapse_and_format(fusion_data), desired


def main():
    parser = argparse.ArgumentParser(
        description="Convert openPBTA fusion table OR list of annofuse files to cbio format."
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
        default="merged_fusion/",
        help="Result output dir. Default is merged_fusion",
    )
    parser.add_argument(
        "-m",
        "--mode",
        action="store",
        dest="mode",
        help="describe source, openX or kfprod, or dgd",
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
    if args.mode != "openX" and args.mode != "kfprod" and args.mode != "dgd":
        print(
            f"-m mode argument must be one of openX, kfprod, or dgd. It is case sensitive. You put {args.mode}",
            file=sys.stderr,
        )
        sys.exit(1)
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    # deal only with RNA metadata
    r_ext: str = "fusion"
    if args.mode == "openX":
        r_ext = "rsem"
    elif args.mode == "dgd":
        r_ext = "DGD_FUSION"
    # ensure sample name is imported as str
    all_file_meta: pd.DataFrame = pd.read_csv(args.table, sep="\t", dtype={"Cbio_Tumor_Name": str})
    # ext used in pbta vs openpedcan varies
    rna_subset: pd.DataFrame = all_file_meta.loc[all_file_meta["File_Type"] == r_ext]
    if rna_subset is None:
        r_ext = "rsem"
        rna_subset = all_file_meta.loc[all_file_meta["File_Type"] == r_ext]
    # reset index so that references work later while iterating
    rna_subset = rna_subset.reset_index(drop=True)
    project_list: np.ndarray = rna_subset.Cbio_project.unique()
    cbio_master, present_cols = init_cbio_master(args.fusion_results, args.mode, rna_subset)

    # Get relevant columns
    cbio_master.set_index("Sample", inplace=True)
    cbio_master.reset_index(inplace=True)
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
    present_rename = {}
    for old, new in rename_dict.items():
        if old in present_cols:
            present_rename[old] = new
    cbio_master.rename(columns=present_rename, inplace=True)
    # Fill in some defaults
    cbio_master["Class"] = "FUSION"
    cbio_master["SV_Status"] = "SOMATIC"
    cbio_master["Site1_Ensembl_Transcript_Id"] = ""
    cbio_master["Site1_Entrez_Gene_Id"] = ""
    cbio_master["Site1_Exon"] = ""
    cbio_master["Site2_Ensembl_Transcript_Id"] = ""
    cbio_master["Site2_Entrez_Gene_Id"] = ""
    cbio_master["Site2_Exon"] = ""
    cbio_master["NCBI_Build"] = "GRCh38"
    cbio_master["Connection_Type"] = "5to3"
    cbio_master["RNA_Support"] = "Yes"
    # Split some columns that have 2 cols worth of info
    if args.mode != "dgd":
        cbio_master[["Site1_Chromosome", "Site1_Position"]] = cbio_master.LeftBreakpoint.str.split(
            ":", expand=True
        )
        cbio_master[["Site2_Chromosome", "Site2_Position"]] = cbio_master.RightBreakpoint.str.split(
            ":", expand=True
        )
    else:
        cbio_master["Site1_Chromosome"] = ""
        cbio_master["Site1_Position"] = ""
        cbio_master["Site2_Chromosome"] = ""
        cbio_master["Site2_Position"] = ""
    # Reformat values to fit needs to be ALL CAPS, replace - with _, remove weird chars
    if args.mode != "dgd":
        cbio_master["Site2_Effect_On_Frame"] = cbio_master["Site2_Effect_On_Frame"].str.upper()
        cbio_master["Site2_Effect_On_Frame"] = cbio_master["Site2_Effect_On_Frame"].str.replace(
            "-", "_"
        )
    else:
        cbio_master["Site2_Effect_On_Frame"] = ""
    # Drop unneeded cols
    cbio_master.drop(["LeftBreakpoint", "RightBreakpoint"], axis=1, inplace=True)

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
    present_order: list[str] = [col for col in order_list if col in cbio_master]
    cbio_master = cbio_master[present_order]
    for project in project_list:
        sub_sample_list = list(
            rna_subset.loc[rna_subset["Cbio_project"] == project, "Cbio_Tumor_Name"]
        )
        fus_fname = out_dir + project + ".fusions.txt"
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
