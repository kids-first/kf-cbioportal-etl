#!/usr/bin/env python3
"""Helper script to append panel data to an existing merged maf file.

Add HGVSP_short if needed and additional fields
"""
import argparse
import csv

import pandas as pd
from add_germline_to_somatic_maf import append_amino_acids_field


def main():

    parser = argparse.ArgumentParser(
        description="Merge panel MAF into existing",
    )
    parser.add_argument(
        "-m", "--main", action="store", help="Main MAF",
    )
    parser.add_argument(
        "-p", "--panel", action="store", help="Panel MAF",
    )
    parser.add_argument(
        "-o",
        "--out_fname",
        action="store",
        help="Result output file name",
    )

    args = parser.parse_args()
    main_maf = pd.read_csv(args.main, sep="\t", skiprows=1, keep_default_na=False, na_values=[""])
    panel_maf = pd.read_csv(args.panel, sep="\t", keep_default_na=False, na_values=[""])
    # immediately drop entrez_id and CNVs
    panel_maf.drop(columns=["Entrez_Gene_Id"], inplace=True)
    panel_maf = panel_maf[(panel_maf["Variant_Type"] != "CNV") & (panel_maf["Tumor_Sample_Barcode"] != "PBCBKE")]

    combined_header = main_maf.columns.to_list()
    combined_header.extend([x for x in panel_maf.columns.to_list() if x not in combined_header])
    # add HGVSP_short if not present
    if "HGVSp_Short" not in panel_maf.columns:
        panel_maf.HGVSp = panel_maf.HGVSp.fillna("")
        panel_maf["HGVSp_Short"] = panel_maf["HGVSp"].apply(append_amino_acids_field)
    merged_maf = pd.concat([main_maf, panel_maf], axis=0)
    # pandas gets cute with making int float, fix that
    int_cols = ["t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count"]
    for col in int_cols:
        if col in merged_maf.columns:
            merged_maf[col] = pd.to_numeric(merged_maf[col], errors="coerce").astype("Int64")
    merged_maf[combined_header].to_csv(args.out_fname, sep="\t", mode="w", index=False, quoting=csv.QUOTE_NONE)


if __name__ == "__main__":

    main()

