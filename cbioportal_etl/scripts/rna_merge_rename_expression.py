#!/usr/bin/env python3
"""Merge rsem gene outputs for cBio usage.

Converts to wide-format
Replace BS ID with cBio ID
Repeat gene names merge to take highest mean expression
Final result also z scored across cohort as log2(FPKM + 1) or as log2(TPM + 1)
"""

import argparse
import concurrent.futures
import io
import json
import os
import sys
import tarfile

import numpy as np
import pandas as pd
from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths
from scipy import stats


def load_rsem_file(rsem_file: str, sample: str, rsem_dir: str, expr_type: str) -> pd.DataFrame:
    """Reads and formats a single RSEM file.
    
    Args: 
        rsem_file: Filename of RSEM (example: 'sample.rsem.genes.results.gz')
        sample: Sample ID
        rsem_dir: Directory where RSEM files are located
        expr_type: Type of expression value to extract (TPM or FPKM)
    """
    try:
        current = pd.read_csv(os.path.join(rsem_dir, rsem_file), sep="\t", index_col=0)
        subset = current[[expr_type]].copy()
        subset.rename(columns={expr_type: sample}, inplace=True)
        return subset
    except Exception as e:
        print(f"{e}\nFailed to process {rsem_file}", file=sys.stderr)
        return None


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
        "-r", 
        "--rsem-dir", 
        action="store", 
        dest="rsem_dir", 
        help="rsem file directory"
    )
    parser.add_argument(
        "-et", 
        "--expression-type", 
        action="store", 
        dest="expression_type", 
        choices=["TPM", "FPKM"], 
        default="TPM", 
        help="Which expression value to use: TPM or FPKM. Default is TPM."
    )
    parser.add_argument(
        "-sc", 
        "--study-config", 
        action="store", 
        dest="study_config", 
        help="cbio study config file."
    )
    parser.add_argument(
        "-dmt", 
        "--default-match-type", 
        action="store", 
        dest="default_match_type", 
        choices=["polyA", "totalRNA", "none"], 
        default="none", 
        help="Default match type for samples with unknown RNA library type for z-score calculations. Use 'polyA' or 'totalRNA' to override fallback to intra-cohort z-score."
    )    
    args = parser.parse_args()

    rsem_dir = args.rsem_dir.rstrip("/")
    out_dir = "merged_rsem/"
    os.makedirs(out_dir, exist_ok=True)

    all_file_meta = pd.read_csv(args.table, sep="\t", dtype={"cbio_sample_name": str})
    rna_subset = all_file_meta[all_file_meta["etl_file_type"] == "rsem"].copy()
    rsem_list = rna_subset[["file_name", "cbio_sample_name"]].drop_duplicates().values.tolist()

    print("Reading RSEM files...", file=sys.stderr)
    seen = set()
    df_list = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(load_rsem_file, fname, sample, rsem_dir, args.expression_type): sample
            for fname, sample in rsem_list if sample not in seen and not seen.add(sample)
        }
        for i, fut in enumerate(concurrent.futures.as_completed(futures), 1):
            result = fut.result()
            if result is not None:
                df_list.append(result)
            if i % 50 == 0:
                print(f"Loaded {i} files", file=sys.stderr)

    master_tbl = pd.concat(df_list, axis=1).astype(np.float32).copy()
    master_tbl.reset_index(inplace=True)
    master_tbl["Hugo_Symbol"] = master_tbl["gene_id"].str.split("_", n=1).str[1]
    master_tbl.drop(columns=["gene_id"], inplace=True)

    sample_list = master_tbl.columns.drop("Hugo_Symbol").tolist()
    master_tbl["mean_expr"] = master_tbl[sample_list].mean(axis=1)
    master_tbl.sort_values("mean_expr", ascending=False, inplace=True)
    master_tbl.drop_duplicates("Hugo_Symbol", keep="first", inplace=True)
    master_tbl.drop(columns=["mean_expr"], inplace=True)
    master_tbl.set_index("Hugo_Symbol", inplace=True)
    log_master_tbl = np.log2(master_tbl + 1)

    project_list = rna_subset.cbio_project.unique()
    print(f"Outputting {args.expression_type} expression results", file=sys.stderr)
    for project in project_list:
        sub_samples = rna_subset[rna_subset["cbio_project"] == project]["cbio_sample_name"].tolist()
        master_tbl[sub_samples].to_csv(f"{out_dir}{project}.rsem_merged.{args.expression_type}.txt", sep="\t")

    print("Calculating z-scores...", file=sys.stderr)
    # Studies with library type column will be processed by library type 
    if "etl_experiment_strategy" in rna_subset.columns:
        # load healthy references
        TOOL_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        with open(args.study_config) as f:
            config_data: dict = json.load(f)
        config_data = resolve_config_paths(config_data, TOOL_DIR)

        healthy_files = {}

        archive_path = config_data["rsem_ref"]["archive"]
        with tarfile.open(archive_path, "r:gz") as tar:
            for member in tar.getmembers():
                if not member.isfile():
                    continue
                name = os.path.basename(member.name)
                if name.endswith(f"_{args.expression_type}.tsv"):
                    key = name.replace(f"_log_{args.expression_type}.tsv", "") + f"_{args.expression_type}"
                    file_obj = tar.extractfile(member)
                    if file_obj:
                        df = pd.read_csv(io.TextIOWrapper(file_obj), sep="\t", index_col=0)
                        # convert genes to Hugo symbols
                        df["Hugo_Symbol"] = df.index.str.split("_", n=1).str[1]
                        # keep the most highly expressed entry per Hugo symbol for duplicates
                        df["mean_expr"] = df.drop(columns="Hugo_Symbol").mean(axis=1)
                        df = df.sort_values("mean_expr", ascending=False).drop_duplicates("Hugo_Symbol", keep="first")
                        df.set_index("Hugo_Symbol", inplace=True)
                        df.drop(columns="mean_expr", inplace=True)
                        healthy_files[key] = df

        zscore_intracohort = []
        zscore_vs_healthy = []
        # Process samples grouped by library type
        unique_strategies = rna_subset["etl_experiment_strategy"].dropna().unique().tolist()
        if rna_subset["etl_experiment_strategy"].isna().any():
            unique_strategies.append(np.nan)
        for library_type in unique_strategies:
            if pd.isna(library_type):
                print("Processing samples with missing etl_experiment_strategy", file=sys.stderr)
                group_df = rna_subset[rna_subset["etl_experiment_strategy"].isna()]
            else:
                print(f"Processing samples with etl_experiment_strategy {library_type}", file=sys.stderr)
                group_df = rna_subset[rna_subset["etl_experiment_strategy"] == library_type]

            group_samples = group_df["cbio_sample_name"].drop_duplicates().tolist()
            if not group_samples:
                continue

            group_tbl = log_master_tbl[group_samples].copy()

            # Intra-cohort z-score
            print(f"Calculating intra-cohort z-score for {library_type}")
            # Calculate z-score and ignore NaN values when computing mean and SD for each row
            # Mean and SD is calculated with non-NaN values
            zscore_vals = stats.zscore(group_tbl, axis=1, nan_policy='omit')
            intracohort_df = pd.DataFrame(zscore_vals, index=group_tbl.index, columns=group_tbl.columns).fillna(0)
            zscore_intracohort.append(intracohort_df)

            # Tumor vs healthy reference z-score (only if match_type available)
            print(f"Calculating tumor vs healthy z-score for {library_type}")
            library_key = str(library_type).strip().lower() if not pd.isna(library_type) else None
            match_type = (
                "totalRNA" if library_key == "rrna depletion"
                else "polyA" if library_key in {"poly-t enrichment", "hybrid selection"}
                else args.default_match_type if args.default_match_type != "none"
                else None
            )
            print(f"  ➤ Match type being used: {match_type}", file=sys.stderr)
            ref_key = f"healthy_{match_type}_{args.expression_type}" if match_type else None

            if ref_key and ref_key in healthy_files:
                healthy_ref = healthy_files[ref_key]
                common_genes = group_tbl.index.intersection(healthy_ref.index)
                mu = healthy_ref.loc[common_genes].mean(axis=1)
                sigma = healthy_ref.loc[common_genes].std(axis=1)
                zscore_df = group_tbl.loc[common_genes].sub(mu, axis=0).div(sigma, axis=0).fillna(0)
                zscore_vs_healthy.append(zscore_df)
            else:
                print(f"No healthy reference available for {library_type} library type, using intra-cohort z-score instead", file=sys.stderr)
                zscore_vs_healthy.append(intracohort_df[group_tbl.columns])

        master_zscore_vs_healthy = pd.concat(zscore_vs_healthy, axis=1).fillna(0).astype(np.float32)
        master_zscore_intracohort = pd.concat(zscore_intracohort, axis=1).fillna(0).astype(np.float32)

        for project in project_list:
            sub_samples = rna_subset[rna_subset["cbio_project"] == project]["cbio_sample_name"].tolist()
            if not master_zscore_vs_healthy.empty:
                healthy_outfile = f"{out_dir}{project}.rsem_merged_vs_healthy_zscore.txt"
                master_zscore_vs_healthy[sub_samples].to_csv(healthy_outfile, sep="\t", float_format="%.4f")

            intra_outfile = f"{out_dir}{project}.rsem_merged_tumor_only_zscore.txt"
            master_zscore_intracohort[sub_samples].to_csv(intra_outfile, sep="\t", float_format="%.4f")

    else:
        # Studies without library type columns will use intra-cohort z-score
        print("No etl_experiment_strategy column found, using intra-cohort z-score", file=sys.stderr)
        zscore_vals = stats.zscore(log_master_tbl, axis=1, nan_policy='omit')
        master_zscore_log = pd.DataFrame(zscore_vals, index=log_master_tbl.index, columns=log_master_tbl.columns).fillna(0).astype(np.float32)

        for project in project_list:
            sub_samples = rna_subset[rna_subset["cbio_project"] == project]["cbio_sample_name"].tolist()
            outfile = f"{out_dir}{project}.rsem_merged_tumor_only_zscore.txt"
            master_zscore_log[sub_samples].to_csv(outfile, sep="\t", float_format="%.4f")
