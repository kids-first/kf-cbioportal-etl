#!/usr/bin/env python3
"""Converts our current fusion and annot SV inputs into the cBio SV format.

For repeat rows (same breakpoint in a sample, different callers),
ARRIBA annot is used with ceiling of mean for counts used
"""

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd
from annotsv_to_cbio_sv import init_cbio_sv_df, setup_outdir_metadata





def main():
    parser = argparse.ArgumentParser(
        description="Convert openPBTA fusion table OR list of annofuse files to cbio format."
    )
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        help="Table with cbio project, kf bs ids, cbio IDs, and file names",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--fusion-results",
        action="store",
        help="annoFuse results dir OR openX merged fusion file",
        required=False,
    )
    parser.add_argument(
        "--sv_results",
        action="store_true",
        help="DNA SV results from annotSV",
        required=False,
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
        help="describe source, openX or kfprod or dgd",
        required=False,
    )
    parser.add_argument(
        "-a",
        "--append",
        action=argparse.BooleanOptionalAction,
        dest="append",
        help="Flag to append, meaning print to STDOUT and skip header",
        required=False,
    )

    args = parser.parse_args()
    print(args.sv_results)
    if args.fusion_results is None and args.sv_results is None:
        print(
            "Either --fusion-results and/or --sv-results must be specified.",
            file=sys.stderr,
        )
        sys.exit(1)

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    # Reorder table
    order_list: list[str] = [
        "Sample_Id",
        "SV_Status",
        "Site1_Hugo_Symbol",
        "Site1_Entrez_Gene_Id",
        "Site1_Ensembl_Transcript_Id",
        "Site1_Exon",
        "Site1_Chromosome",
        "Site1_Contig",
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
        "Tumor_Paired_End_Read_Count",
        "Annotation",
        "DNA_Support",
        "RNA_Support",
        "SV_Length",
        "Connection_Type",
        "Breakpoint_Type",
        "Event_Info",
        "Class",
        "External_Annotation",
        "Comments",
    ]
    if args.sv_results:
        link_input_dir = "annotSV_results"
        dna_sv_subset: pd.DataFrame = setup_outdir_metadata(link_input_dir, args.out_dir, args.table)
        project_list: np.ndarray = dna_sv_subset.cbio_project.unique()
        cbio_sv_df: pd.DataFrame = init_cbio_sv_df(link_input_dir, dna_sv_subset)


        for project in project_list:
            cbio_sv_fname = os.path.join(args.out_dir, project + ".svs.txt")
            cbio_sv_df.set_index("Sample_Id", inplace=True)
            cbio_sv_df.to_csv(cbio_sv_fname, sep="\t", mode="w", index=True, quoting=csv.QUOTE_NONE)

if __name__ == "__main__":
    main()
