#!/usr/bin/env python
"""Convert protein modification data to cBio format.

Select samples only in manifest, reformat columns and rename samples aacordingly, output cBio formatted raw and z scored results
"""
import argparse
import csv
import os

import numpy as np
import pandas as pd
from scipy import stats
from format_total_proteome import get_metadata


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
        "-p",
        "--protein-file",
        action="store",
        help="Protein modification file",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--protein-mod",
        action="store",
        help="Protein modification modification, i.e. phoshp, acetyl, etc",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        action="store",
        dest="out_dir",
        default="merged_prot",
        help="Result output dir. Default is merged_prot",
    )

    args = parser.parse_args()
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    prot_subset: pd.DataFrame = get_metadata(args.table, "EXT_PROT")
    project_list: np.ndarray = prot_subset.cbio_project.unique()
    prot_data: pd.DataFrame = pd.read_csv(
        args.protein_file, sep="\t", keep_default_na=False, na_values=[""]
    )
    rename_index: dict[str, str] = {"id": "ENTITY_STABLE_ID", "id.short": "NAME", "gene_transcript_protein_ENSEMBL": "DESCRIPTION"}
    prot_data.rename(columns=rename_index, inplace=True)
    cbio_sample_name: list[str] = prot_subset["cbio_sample_name"].to_list()
    # until we have BS IDs, will use collection event IDs splitting cBio IDs
    collect_event_id: list[str] = [x.split("_")[0] for x in cbio_sample_name]
    # rename as cBio
    rename_cols: dict[str, str] = dict(zip(collect_event_id, cbio_sample_name, strict=True))
    prot_data.rename(columns=rename_cols, inplace=True)
    prot_data.set_index(["ENTITY_STABLE_ID", "NAME", "DESCRIPTION"], inplace=True)
    # take log2
    prot_data = np.log2(prot_data[cbio_sample_name] + 1).round(4)
    out_fname: str = os.path.join(out_dir, f"{project_list[0]}_{args.protein_mod}_protein.txt")
    prot_data.to_csv(
                out_fname, sep="\t", index=True, quoting=csv.QUOTE_NONE, header=True
            )
    # output z scores
    zscore_vals = stats.zscore(prot_data, axis=1, nan_policy='omit')
    prot_zscore = pd.DataFrame(zscore_vals, index=prot_data.index, columns=prot_data.columns).fillna(0).round(4)
    zscore_fname: str = os.path.join(out_dir, f"{project_list[0]}_{args.protein_mod}_protein_zscore.txt")
    prot_zscore.to_csv(
        zscore_fname, sep="\t", index=True, quoting=csv.QUOTE_NONE, header=True
    )
if __name__ == "__main__":
    main()
