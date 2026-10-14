#!/usr/bin/env python
"""Convert raw total proteome data to cBio format.

Select samples only in manifest, reformat columns and rename samples aacordingly, output cBio formatted raw and z scored results
"""
import argparse
import csv
import os

import numpy as np
import pandas as pd
from scipy import stats


def get_metadata(table: str, ext: str) -> pd.DataFrame:
    """Get metadata for the given file type.

    Args:
        table (str): Path to the metadata table.
        ext (str): File type to filter on.

    Returns:
        pd.DataFrame: Filtered metadata DataFrame.

    """
    # ensure sample name is imported as str
    all_file_meta: pd.DataFrame = pd.read_csv(table, sep="\t", dtype={"cbio_sample_name": str})
    prot_subset: pd.DataFrame = all_file_meta.loc[all_file_meta["etl_file_type"] == ext]
    return prot_subset


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
        "--total-proteome",
        action="store",
        help="Raw proteome file",
        required=False,
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        action="store",
        dest="out_dir",
        default="merged_prot",
        help="Result output dir. Default is merged_prot",
    )
    parser.add_argument(
        "-l",
        "--log2-transform",
        action="store_true",
        help="Flag to log2 transform the data",
        required=False,
    )


    args = parser.parse_args()
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    prot_subset: pd.DataFrame = get_metadata(args.table, "EXT_PROT")
    project_list: np.ndarray = prot_subset.cbio_project.unique()
    prot_data: pd.DataFrame = pd.read_csv(
        args.total_proteome, sep="\t", keep_default_na=False, na_values=[""]
    )
    prot_data['Composite.Element.REF'] = prot_data['gene_name_V42'] + " " + prot_data['gene_transcript_protein_ENSEMBL']
    cbio_sample_name: list[str] = prot_subset["cbio_sample_name"].to_list()
    # until we have BS IDs, will use collection event IDs splitting cBio IDs
    collect_event_id: list[str] = [x.split("_")[0] for x in cbio_sample_name]
    # rename as cBio
    rename_cols: dict[str, str] = dict(zip(collect_event_id, cbio_sample_name, strict=True))
    prot_data.rename(columns=rename_cols, inplace=True)
    prot_data.set_index("Composite.Element.REF", inplace=True)
    # take log2
    if args.log2_transform:
        prot_data = np.log2(prot_data[cbio_sample_name] + 1).round(4)
    else:
        prot_data = prot_data[cbio_sample_name].round(4)
    out_fname: str = os.path.join(out_dir, project_list[0] +"_total_protein.txt")
    prot_data.to_csv(
                out_fname, sep="\t", index=True, quoting=csv.QUOTE_NONE, header=True
            )
    # output z scores
    zscore_vals = stats.zscore(prot_data, axis=1, nan_policy='omit')
    prot_zscore = pd.DataFrame(zscore_vals, index=prot_data.index, columns=prot_data.columns).fillna(0).round(4)
    zscore_fname: str = os.path.join(out_dir, project_list[0] +"_total_protein_zscore.txt")
    prot_zscore.to_csv(
        zscore_fname, sep="\t", index=True, quoting=csv.QUOTE_NONE, header=True
    )
if __name__ == "__main__":
    main()
