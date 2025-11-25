#!/usr/bin/env python3
"""Collapse counts table and Z score cBio usage.

Repeat gene names merge to take highest mean expression
Final result also z scored across cohort as log2(COUNTS + 1)
"""

import argparse
import sys

import numpy as np
import pandas as pd
from scipy import stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Custom collapse coiunts table and calculate z score.")
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        dest="table",
        help="Counts table",
    )
  
    args = parser.parse_args()

    master_tbl = pd.read_csv(args.table, sep="\t", header=0, index_col=0)
    master_tbl.reset_index(inplace=True)
    sample_list = master_tbl.columns.drop("Hugo_Symbol").tolist()
    master_tbl["mean_expr"] = master_tbl[sample_list].mean(axis=1)
    master_tbl.sort_values("mean_expr", ascending=False, inplace=True)
    master_tbl.drop_duplicates("Hugo_Symbol", keep="first", inplace=True)
    master_tbl.drop(columns=["mean_expr"], inplace=True)
    master_tbl.set_index("Hugo_Symbol", inplace=True)
    log_master_tbl = np.log2(master_tbl + 1)

    print(f"Outputting counts expression results", file=sys.stderr)
    master_tbl.to_csv(f"counts_collapsed.txt", sep="\t")

    print("Calculating z-scores...", file=sys.stderr)

    # Studies without library type columns will use intra-cohort z-score
    print("No etl_experiment_strategy column found, using intra-cohort z-score", file=sys.stderr)
    zscore_vals = stats.zscore(log_master_tbl, axis=1, nan_policy='omit')
    master_zscore_log = pd.DataFrame(zscore_vals, index=log_master_tbl.index, columns=log_master_tbl.columns).fillna(0).astype(np.float32)

    master_zscore_log.to_csv("counts_collapsed_zscored.txt", sep="\t", float_format="%.4f")
