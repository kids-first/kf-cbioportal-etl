#!/usr/bin/env python3

import sys
import argparse
import os
import pandas as pd
import numpy as np
import pdb


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
            # concat annofuse file, rename Sample Column according to cBio name
            pdb.set_trace()
            ann_file = pd.read_csv(fusion_results + "/" + flist[i], sep="\t")
            metadata_row = rna_metadata.loc[rna_metadata["File_Name"] == flist[i]]
            ann_file = ann_file.assign(Samplle=metadata_row["Sample"])
            frame_list.append(ann_file)
        concat_frame = pd.concat(frame_list)
        del frame_list
        fusion_data = concat_frame[
            ["Gene1A", "Sample", "FusionName", "Caller", "Fusion_Type"]
        ]
        del concat_frame
        fusion_data = (
            fusion_data.groupby(["Gene1A", "Sample", "FusionName", "Fusion_Type"])[
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
    rna_subset.rename(columns={"T_CL_BS_ID": "Sample"}, inplace=True)
    rna_subset.set_index("Sample", inplace=True)
    project_list = rna_subset.Cbio_project.unique()
    cbio_master = init_cbio_master(args.fusion_results, args.mode, rna_subset)

    # Get relevant columns
    cbio_master.set_index("Sample", inplace=True)
    # get cbio IDs and seq center
    cbio_master = pd.merge(
        cbio_master, rna_subset[["Cbio_Tumor_Name"]], on="Sample", how="inner"
    )
    # drop bs ids
    cbio_master.reset_index(inplace=True)
    cbio_master.drop("Sample", axis=1, inplace=True)

    cbio_master.rename(
        columns={
            "Gene1A": "Hugo_Symbol",
            "FusionName": "Fusion",
            "Caller": "Method",
            "Fusion_Type": "Frame",
            "Cbio_Tumor_Name": "Tumor_Sample_Barcode",
            "SQ_Value": "Center",
        },
        inplace=True,
    )

    cbio_master["DNA_support"] = "no"
    cbio_master["RNA_support"] = "yes"
    # create blank entrez ID column so that 3' gene names can be searched
    cbio_master["Entrez_Gene_Id"] = ""
    # Reorder table
    order_list = [
        "Sample_ID",
        "Site1_Hugo_Symbol",
        "Site1_Entrez_Gene_Id",
        "Site1_Chromosome",
        "Site1_Position",
        "Site2_Hugo_Symbol",
        "Site2_Entrez_Gene_Id",
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
    cbio_master.set_index("Hugo_Symbol", inplace=True)
    for project in project_list:
        sub_sample_list = list(
            rna_subset.loc[rna_subset["Cbio_project"] == project, "Cbio_Tumor_Name"]
        )
        fus_fname = out_dir + project + ".fusions.txt"
        fus_file = open(fus_fname, "w")
        fus_tbl = cbio_master[cbio_master.Tumor_Sample_Barcode.isin(sub_sample_list)]
        fus_tbl.reset_index(inplace=True)
        fus_tbl.fillna("NA", inplace=True)
        fus_head = "\t".join(fus_tbl.columns) + "\n"
        fus_file.write(fus_head)
        fus_data = list(fus_tbl.values)
        # "hack" to allow 3' end of fusion to be searched
        for data in fus_data:
            fus_file.write("\t".join(data) + "\n")
            (a, b) = data[4].split("--")
            if a != b:
                data[0] = b
                fus_file.write("\t".join(data) + "\n")
        fus_file.close()
