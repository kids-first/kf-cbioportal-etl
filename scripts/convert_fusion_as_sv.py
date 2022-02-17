#!/usr/bin/env python3
"""
This is a beta script for converting our current fusion inputs into the new cBio SV format.
However, the front end support is not complete, so this will not be used until it is.
"""


import sys
import argparse
import os
import pandas as pd
import pyranges
import csv


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert openPBTA fusion table OR list of annofuse files to cbio format."
    )
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        dest="table",
        help="Table with cbio project, kf bs ids, cbio IDs, and file names",
    )
    parser.add_argument(
        "-f",
        "--fusion-results",
        action="store",
        dest="fusion_results",
        help="annoFuse results dir",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        action="store",
        dest="out_dir",
        default="merged_fusion/",
        help="Result output dir. Default is merged_fusion",
    )
    args = parser.parse_args()

    def init_cbio_master(fusion_results, rna_metadata):
        """
        Use data frame subset on RNA fusion files to find and merge result files
        """
        flist = rna_metadata.File_Name
        frame_list = []
        for i in range(0, len(flist), 1):
            try:
                # concat annofuse file, rename Sample Column according to cBio name
                ann_file = pd.read_csv(fusion_results + "/" + flist[i], sep="\t")
                ann_file = ann_file.assign(Sample=rna_metadata.at[i, "Cbio_Tumor_Name"])
                frame_list.append(ann_file)
            except Exception as e:
                pdb.set_trace()
                hold = 1
        concat_frame = pd.concat(frame_list)
        del frame_list
        fusion_data = concat_frame[
            [
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
                "Caller"
                ]
        ]
        del concat_frame
        fusion_data = (
            fusion_data.groupby([
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
                "Fusion_anno",])[
                "Caller"
            ]
            .apply(", ".join)
            .reset_index()
        )
        fusion_data["Caller"] = fusion_data["Caller"].str.upper()
        return fusion_data

    out_dir = args.out_dir
    try:
        os.mkdir(out_dir)
    except:
        sys.stderr.write("output dir already exists\n")
        sys.stderr.flush()

    # deal only with RNA metadata
    r_ext = "fusion"
    all_file_meta = pd.read_csv(args.table, sep="\t")
        
    rna_subset = all_file_meta.loc[all_file_meta["File_Type"] == r_ext]
    # reset index so that references work later while iterating
    rna_subset = rna_subset.reset_index(drop=True)
    project_list = rna_subset.Cbio_project.unique()
    cbio_master = init_cbio_master(args.fusion_results, rna_subset)

    # Get relevant columns
    cbio_master.set_index("Sample", inplace=True)
    # get cbio IDs and seq center

    # drop bs ids
    cbio_master.reset_index(inplace=True)

    cbio_master.rename(
        columns={
            "Sample": "Sample_ID",
            "Gene1A": "Site1_Hugo_Symbol",
            "Gene1B": "Site2_Hugo_Symbol",
            "Fusion_Type": "Site2_Effect_On_Frame",
            "JunctionReadCount": "Tumor_Read_Count",
            "SpanningFragCount": "Tumor_Split_Read_Count",
            "annots": "Annotation",
            "FusionName": "Event_Info",
            "Fusion_anno": "External_Annotation",
            "Caller": "Comments"
        },
        inplace=True
    )
    # Fill in some defaults
    cbio_master['Class'] = "FUSION"
    cbio_master["Site1_Ensembl_Transcript_Id"] = ""
    cbio_master["Site1_Entrez_Gene_Id"] = ""
    cbio_master["Site1_Exon"] = ""
    cbio_master["Site2_Ensembl_Transcript_Id"] = ""
    cbio_master["Site2_Entrez_Gene_Id"] = ""
    cbio_master["Site2_Exon"] = ""
    cbio_master["NCBI_Build"] = "GRCh38"
    cbio_master["Connection_Type"] = "5to3"
    # Split some columns that have 2 cols worth of info
    cbio_master[['Site1_Chromosome','Site1_Position']] = cbio_master.LeftBreakpoint.str.split(":", expand=True)
    cbio_master[['Site2_Chromosome','Site2_Position']] = cbio_master.RightBreakpoint.str.split(":", expand=True)
    # cbio_master['Class'] = cbio_master.Annotation.str.split("],").str[1]
    # cbio_master['Annotation'] = cbio_master.Annotation.str.split("],").str[0]
    # Reformat values to fit needs to be ALL CAPS, replace - with _, remove weird chars
    cbio_master["Site2_Effect_On_Frame"] = cbio_master["Site2_Effect_On_Frame"].str.upper()
    cbio_master['Site2_Effect_On_Frame'] = cbio_master['Site2_Effect_On_Frame'].str.replace('-','_')
    # cbio_master.loc[cbio_master["Annotation"] == "[", "Annotation"] = "NA"
    # cbio_master['Class'] = cbio_master['Class'].str.upper()
    # Drop unneeded cols
    cbio_master.drop(['LeftBreakpoint', 'RightBreakpoint'], axis=1, inplace=True)

    # Reorder table
    order_list = [
        "Sample_ID",
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
        "Connection_Type",
        "Event_Info",
        "Class",
        "External_Annotation"
    ]
    cbio_master = cbio_master[order_list]
    # cbio_master.set_index("Sample_ID", inplace=True)
    for project in project_list:
        sub_sample_list = list(
            rna_subset.loc[rna_subset["Cbio_project"] == project, "Cbio_Tumor_Name"]
        )
        fus_fname = out_dir + project + ".fusions.txt"
        fus_tbl = cbio_master[cbio_master.Sample_ID.isin(sub_sample_list)]
        fus_tbl.fillna("NA", inplace=True)
        fus_tbl.set_index("Sample_ID", inplace=True)
        # fus_head = "\t".join(fus_tbl.columns) + "\n"
        # fus_file.write(fus_head)
        # fus_data = list(fus_tbl.values)
        # # "hack" to allow 3' end of fusion to be searched
        # for data in fus_data:
        #     fus_file.write("\t".join(data) + "\n")
        #     (a, b) = data[4].split("--")
        #     if a != b:
        #         data[0] = b
        #         fus_file.write("\t".join(data) + "\n")
        fus_tbl.to_csv(fus_fname, sep="\t", mode="w", index=True, quoting=csv.QUOTE_NONE)
