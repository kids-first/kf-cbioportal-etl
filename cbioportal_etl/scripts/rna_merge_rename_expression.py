#!/usr/bin/env python3
"""Merge rsem gene outputs for cBio usage.

Converts to wide-format
Replace BS ID with cBio ID
Repeat gene names merge to take highest mean expression
Final result also z scored across cohort as log2(FPKM + 1) or as log2(TPM + 1)
"""

import argparse
import concurrent.futures
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats


def mt_collate_df(rsem_file: str) -> int:
    """Concat desired FPKM/TPM data into data frame.

    Args:
        rsem_file: File name of current RSEM file to process
    """
    try:
        sample: str = rna_subset.loc[rna_subset["file_name"] == rsem_file, "cbio_sample_name"].iloc[
            0
        ]
        if sample not in seen_dict:
            current: pd.DataFrame = pd.read_csv(rsem_dir + rsem_file, sep="\t", index_col=0)
            cur_subset = current[[args.expression_type]].copy()
            cur_subset.rename(columns={args.expression_type: sample}, inplace=True)
            df_list.append(cur_subset)
            seen_dict[sample] = 1
        else:
            print(f"{sample} was already seen, likely due to multiple dx, skipping!")

    except Exception as e:
        print(f"{e}\nFailed to process {rsem_file}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge rsem files using cavatica file info.")
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        dest="table",
        help="Table with cbio project, kf bs ids, cbio IDs, and file names",
    )
    parser.add_argument(
        "-r", "--rsem-dir", action="store", dest="rsem_dir", help="rsem file directory"
    )
    parser.add_argument(
        "-ds", "--data-sample", action="store", dest="data_sample", default="datasheets/data_clinical_sample.txt", help="data_clinical_sample.txt file"
    )
    parser.add_argument(
        "-et", "--expression-type", action="store", dest="expression_type", choices=["TPM", "FPKM"], default="FPKM", help="Which expression value to use: TPM or FPKM. Default is FPKM."
    )
    args = parser.parse_args()

    rsem_dir: str = args.rsem_dir

    if rsem_dir[-1] != "/":
        rsem_dir += "/"

    # file_meta_dict = get_file_metadata(args.table, 'rsem')
    out_dir = "merged_rsem/"
    os.makedirs(out_dir, exist_ok=True)

    # Use cbio table to create master table to choose representative gene symbol, calc z scores, output both
    # ensure sample name is imported as str
    all_file_meta: pd.DataFrame = pd.read_csv(args.table, sep="\t", dtype={"cbio_sample_name": str})
    rna_subset: pd.DataFrame = all_file_meta.loc[all_file_meta["etl_file_type"] == "rsem"]
    rsem_list: list[str] = rna_subset["file_name"].to_list()
    print("Creating merged rsem table", file=sys.stderr)
    x: int = 1
    m: int = 50

    df_list: list[pd.DataFrame] = []
    # some samples have more than one dx, need only process once
    seen_dict: dict[str, int] = {}
    with concurrent.futures.ThreadPoolExecutor(16) as executor:
        results = {
            executor.submit(mt_collate_df, rsem_list[i]): rsem_list[i]
            for i in range(len(rsem_list))
        }
        for result in concurrent.futures.as_completed(results):
            if x % m == 0:
                print(f"Merged {x} files", file=sys.stderr)
                sys.stderr.flush()
            x += 1
    master_tbl = pd.concat(df_list, axis=1)
    del df_list

    print("Merge complete, picking representative gene transcripts", file=sys.stderr)
    sys.stderr.flush()
    master_tbl.reset_index(inplace=True)
    # drop ens_id_sym and replace with just sym
    gene_id_sym_list: list[str] = master_tbl["gene_id"].to_list()
    gene_sym_list = []
    for entry in gene_id_sym_list:
        parts: list[str] = entry.split("_")
        gene_sym = "_".join(parts[1:])
        gene_sym_list.append(gene_sym)
    master_tbl["Hugo_Symbol"] = gene_sym_list
    master_tbl.drop(columns=["gene_id"], inplace=True)
    sample_list: list[str] = list(master_tbl.columns.values)
    # remove value with Hugo_Symbol
    h_idx: int = sample_list.index("Hugo_Symbol")
    sample_list.pop(h_idx)
    # look for repeated instances of gene symbol, use highest average exp as rep value
    dup_hugo_tbl: pd.DataFrame = master_tbl[master_tbl.duplicated(["Hugo_Symbol"])]
    rpt_sym: np.ndarray = dup_hugo_tbl.Hugo_Symbol.unique()
    for sym in rpt_sym:
        gene_eval: pd.DataFrame = master_tbl.loc[master_tbl["Hugo_Symbol"] == sym]
        mean_expr: pd.Series = gene_eval[sample_list].mean(axis=1).sort_values(ascending=False)
        to_dump = list(mean_expr.index)[1:]
        master_tbl.drop(to_dump, inplace=True)
    master_tbl.set_index("Hugo_Symbol", inplace=True)
    gene_sym_list = master_tbl.index
    project_list = rna_subset.cbio_project.unique()

    print(f"Outputting {args.expression_type} expression results", file=sys.stderr)
    sys.stderr.flush()
    for project in project_list:
        sub_sample_list = list(
            rna_subset.loc[rna_subset["cbio_project"] == project, "cbio_sample_name"]
        )
        expr_fname = f"{out_dir}{project}.rsem_merged.{args.expression_type}.txt"
        master_tbl[sub_sample_list].to_csv(expr_fname, sep="\t", mode="w", index=True)

    print("Calculating z scores", file=sys.stderr)
    sys.stderr.flush()

    clinical_df = pd.read_csv(args.data_sample, sep="\t", comment="#")
    if "RNA_LIBRARY_SELECTION" in clinical_df.columns:
        print("Detected RNA_LIBRARY_SELECTION column, calculating group-specific z-scores", file=sys.stderr)
        clinical_df = clinical_df[["SAMPLE_ID", "RNA_LIBRARY_SELECTION"]].drop_duplicates()
        rna_subset = rna_subset.merge(clinical_df, left_on="cbio_sample_name", right_on="SAMPLE_ID", how="left")

        healthy_ref_dir = "/home/ubuntu/tools/kf-cbioportal-etl/cbioportal_etl/REFS"
        healthy_files = {
            "polyA_FPKM": pd.read_csv(f"{healthy_ref_dir}/healthy_polyA_log_FPKM.tsv", sep="\t", index_col=0),
            "polyA_TPM": pd.read_csv(f"{healthy_ref_dir}/healthy_polyA_log_TPM.tsv", sep="\t", index_col=0),
            "totalRNA_FPKM": pd.read_csv(f"{healthy_ref_dir}/healthy_totalRNA_log_FPKM.tsv", sep="\t", index_col=0),
            "totalRNA_TPM": pd.read_csv(f"{healthy_ref_dir}/healthy_totalRNA_log_TPM.tsv", sep="\t", index_col=0),
        }

        group_zscores = []
        group_zscores_vs_healthy = []

        for library_type in rna_subset.RNA_LIBRARY_SELECTION.unique():
            print(f"Calculating z-scores for RNA library type: {library_type}", file=sys.stderr)
            if pd.isna(library_type):
                group_samples = rna_subset[rna_subset["RNA_LIBRARY_SELECTION"].isna()]["cbio_sample_name"].tolist()
            else:
                group_samples = rna_subset[rna_subset["RNA_LIBRARY_SELECTION"] == library_type]["cbio_sample_name"].tolist()

            if not group_samples:
                continue

            group_tbl = master_tbl[group_samples]

            # Intra-cohort z-score
            group_z = stats.zscore(np.log2(np.array(group_tbl + 1)), axis=1)
            group_z_df = pd.DataFrame(group_z, index=master_tbl.index, columns=group_samples)
            group_zscores.append(group_z_df)

            # Tumor vs healthy z-score
            lib_key = str(library_type).strip().lower()
            if lib_key == "rrna depletion":
                match_type = "totalRNA"
            elif lib_key in ["poly-t enrichment", "hybrid selection"]:
                match_type = "polyA"
            else:
                match_type = None
        
            # Construct key based on expression type
            if match_type:
                ref_key = f"{match_type}_{args.expression_type}"
                if ref_key in healthy_files:
                    healthy_ref = healthy_files[ref_key]
                    healthy_mu = healthy_ref.mean(axis=1)
                    healthy_sigma = healthy_ref.std(axis=1)

                    group_tbl_log = np.log2(group_tbl + 1)
                    z_vs_healthy = group_tbl_log.sub(healthy_mu, axis=0).div(healthy_sigma, axis=0)
                    group_zscores_vs_healthy.append(z_vs_healthy)
                else:
                    print(f"No matching healthy reference for library type {library_type}, calculating intra-cohort z-score", file=sys.stderr)
                    z_fallback = stats.zscore(np.log2(np.array(group_tbl + 1)), axis=1)
                    group_zscores_vs_healthy.append(pd.DataFrame(z_fallback, index=master_tbl.index, columns=group_samples))
            else:
                print(f"No library type '{library_type}', calculating intra-cohort z-score", file=sys.stderr)
                z_fallback = stats.zscore(np.log2(np.array(group_tbl + 1)), axis=1)
                group_zscores_vs_healthy.append(pd.DataFrame(z_fallback, index=master_tbl.index, columns=group_samples))

        # concat the intra-cohort zscores
        master_zscore_log = pd.concat(group_zscores, axis=1)
        master_zscore_log.fillna(0, inplace=True)
        # concat the tumor vs healthy zscores and output it
        master_zscore_vs_healthy = pd.concat(group_zscores_vs_healthy, axis=1)
        master_zscore_vs_healthy.fillna(0, inplace=True)
        for project in project_list:
            sub_sample_list = list(
                rna_subset.loc[rna_subset["cbio_project"] == project, "cbio_sample_name"]
            )
            zscore_healthy_fname = f"{out_dir}{project}.rsem_merged_vs_healthy_zscore.txt"
            master_zscore_vs_healthy[sub_sample_list].to_csv(
                zscore_healthy_fname,
                sep="\t",
                mode="w",
                index=True,
                float_format="%.4f",
            )
    else: 
        print("No RNA_LIBRARY_SELECTION column found, using full cohort for z-score calculation", file=sys.stderr)
        z_scored = stats.zscore(np.log2(np.array(master_tbl + 1)), axis=1)
        del master_tbl
        master_zscore_log: pd.DataFrame = pd.DataFrame(
            z_scored, index=gene_sym_list, columns=sample_list
        )
    # this may be memory-intensive for some insane reason...
    print("Replacing NaN with 0", file=sys.stderr)
    sys.stderr.flush()
    master_zscore_log.fillna(0, inplace=True)
    print("Outputting z scored results", file=sys.stderr)
    sys.stderr.flush()

    for project in project_list:
        sub_sample_list = list(
            rna_subset.loc[rna_subset["cbio_project"] == project, "cbio_sample_name"]
        )
        zscore_fname = f"{out_dir}{project}.rsem_merged_tumor_only_zscore.txt"
        master_zscore_log[sub_sample_list].to_csv(
            zscore_fname,
            sep="\t",
            mode="w",
            index=True,
            float_format="%.4f",
        )
