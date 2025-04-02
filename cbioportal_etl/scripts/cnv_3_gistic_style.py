#!/usr/bin/env python3
"""Convert wide-format CNV values to discrete GISTIC-style wide-format table."""

import argparse
import json
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from get_file_metadata_helper import get_file_metadata

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths

parser = argparse.ArgumentParser(
    description="Convert merged cnv values to discrete coded values.",
)
parser.add_argument(
    "-d",
    "--merged_cnv_dir",
    action="store",
    dest="merged_cnv_dir",
    help="merged cnv dir",
)
parser.add_argument(
    "-j",
    "--config",
    action="store",
    dest="config_file",
    help="json config file with data types and data locations",
)
parser.add_argument(
    "-i",
    "--info_dir",
    action="store",
    dest="info_dir",
    help="If info files available, provide here",
)
parser.add_argument(
    "-t",
    "--table",
    action="store",
    dest="table",
    help="Table with cbio project, kf bs ids, cbio IDs, and file names",
)


def adjust_cn(bs_id: str) -> int:
    """Adjust GISTIC Value based on sample ploidy.

    Need to use ploidy info to ensure as sample with ploidy 3 for instance,
    has CN value of 3 be 0 and adjust all other accordingly
    Args:
        bs_id: BS ID of sample currently being considered
    Returns: list of int representing success or failure and BS ID
    """
    try:
        # If there is a info file manifest, get ploidy. Else skip and assume ploidy 2
        samp_id: str = bs_cbio_dict[bs_id]
        ploidy: int = 2

        if manifest is not None:
            with open(args.info_dir + "/" + manifest.loc[bs_id]["File_Name"]) as info:
                for datum in info:
                    datum = datum.rstrip("\n")
                    try:
                        if datum.startswith("Output_Ploidy"):
                            (key, value) = datum.split("\t")
                            ploidy = int(value)
                            break
                    except Exception as e:
                        print(f"WARN: {e} entry could not be split", file=sys.stderr)
        data[samp_id] = data[samp_id] - ploidy
        min_cn: int = ploidy * -1
        # adjust discrete low range only in ploidy > 2
        if min_cn != DEFAULT_MIN_CN:
            data[samp_id] = np.where(data[samp_id].between((min_cn + 1), -2), -1, data[samp_id])
            data[samp_id] = np.where(data[samp_id] == min_cn, -2, data[samp_id])
        data[samp_id] = np.where(data[samp_id].between(2, high_gain), 1, data[samp_id])
        data[samp_id] = np.where(data[samp_id] > high_gain, 2, data[samp_id])
        return 0

    except Exception as E:
        print(f"{E}", file=sys.stderr)
        print(f"Error processing sample entry for {bs_id}\n", file=sys.stderr)
        return 1


args = parser.parse_args()
TOOL_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

with open(args.config_file) as f:
    config_data: dict = json.load(f)
config_data = resolve_config_paths(config_data, TOOL_DIR)


flist: bytes = subprocess.check_output(
    "find " + args.merged_cnv_dir + " -name *.predicted_cnv.txt", shell=True
)

fname_list: list[str] = flist.decode().split("\n")
file_meta_dict: dict[str, dict[str, dict[str, str]]] = get_file_metadata(args.table, "cnv")

manifest = None
if args.info_dir is not None:
    manifest: pd.DataFrame = pd.read_csv(args.table, sep="\t")
    manifest.set_index(["T_CL_BS_ID"], inplace=True)
    manifest = manifest.loc[manifest["File_Type"] == "info"]

else:
    print("No info file given, will assume ploidy 2\n", file=sys.stderr)

if fname_list[-1] == "":
    fname_list.pop()
for fname in fname_list:
    parts: Match[str] | None = re.search(
        "^" + args.merged_cnv_dir + "/(.*).predicted_cnv.txt", fname
    )
    cbio_dx: str = parts.group(1)
    try:
        data: pd.DataFrame = pd.read_csv(fname, sep="\t")
    except Exception as e:
        print(e, file=sys.stderr)
    data.set_index("Hugo_Symbol")
    # sample list would be cbio ids
    samp_list: list[str] = list(data.columns)[1:]
    bs_cbio_dict: dict[str, str] = {}
    for samp_id in samp_list:
        bs_id: str = file_meta_dict[cbio_dx][samp_id]["kf_tum_id"]
        bs_cbio_dict[bs_id] = samp_id
    high_gain: int = config_data["cnv_high_gain"]

    x: int = 1
    m: int = 50
    DEFAULT_MIN_CN = -2
    for bs_id in bs_cbio_dict:
        exit_code = adjust_cn(bs_id)
        if exit_code == 1:
            print(f"Had trouble processing data for {bs_id}", file=sys.stderr)
            sys.exit(1)
        if x % m == 0:
            print("Processed {x} samples", file=sys.stderr)
            sys.stderr.flush()
        x += 1

    print("Conversion completed. Writing results to file", file=sys.stderr)
    new_fname = f"{args.merged_cnv_dir}/{cbio_dx}.discrete_cnvs.txt"
    data.to_csv(new_fname, sep="\t", mode="w", index=False)
