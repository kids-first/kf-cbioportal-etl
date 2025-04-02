#!/usr/bin/env python3
"""Merges files in gene <tab> entrez id <tab> copy number format into a genes-by-sample copy number table.

Uses info files to get ploidy (set neutral value) and ETL cbio manifest to get file names
"""

import argparse
import concurrent.futures
import json
import os
import sys
from typing import IO

import pandas as pd
from get_file_metadata_helper import get_file_metadata

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths


def process_cnv(cnv_fn: str, cur_cnv_dict: dict, samp_id: str) -> dict[str, dict[str, str]]:
    """Initialize cnv dict with gene and samp ID as keys, CNV as value.

    Args:
        cnv_fn: string of CNV file name to process
        cur_cnv_dict: Dict to house CNV data
        samp_id: cBio sample ID associated with this file
    Returns:
        Dict with gene and sample ID as key, copy number as value
    """
    for entry in open(cnv_fn):
        (gene, _gid, value) = entry.rstrip("\n").split("\t")
        if gene not in cur_cnv_dict:
            cur_cnv_dict[gene] = {}
        if samp_id in cur_cnv_dict[gene]:
            print(
                f"ERROR: Sample ID {samp_id} already processed.  Back to the drawing board!",
                file=sys.stderr,
            )
            exit(1)
        else:
            cur_cnv_dict[gene][samp_id] = value
    return cur_cnv_dict


def get_ploidy(cfree_info_fname: str) -> str:
    """Process ControFreeC Info file for [loidy information.

    Used to set a default value for CN if sample does not have a one for a specific gene
    Args:
        cfree_info_fname: Filename of ControFreeC info file
    Returns: str representation of output ploidy
    """
    info: IO = open(cfree_info_fname)
    for datum in info:
        datum: str = datum.rstrip("\n")
        try:
            if datum.startswith("Output_Ploidy"):
                (_key, value) = datum.split("\t")
                return value
        except Exception as e:
            print(f"WARN: {e} entry could not be split", file=sys.stderr)
    return "2"  # default if ploidy can't be found


def process_table(
    cbio_dx: str, file_meta_dict: dict[str, dict[str, dict[str, str]]]
) -> tuple[int, str]:
    """Create CNV table in wide format.

    Args:
        cbio_dx: Disease/project name to process
        file_meta_dict: Dict from file metadata
    Returns:
        tuple of int representing success or failure, cBio Sample ID
    """
    try:
        # project/disease name should be name of directory hosting datasheet
        print("Processing " + cbio_dx + " project", file=sys.stderr)
        new_cnv: IO = open(out_dir + cbio_dx + ".predicted_cnv.txt", "w")
        cur_cnv_dict: dict = {}
        s_list: list = []
        ploidy_dict: dict = {}
        for cbio_tum_id in file_meta_dict[cbio_dx]:
            orig_fname: str = file_meta_dict[cbio_dx][cbio_tum_id]["fname"]
            kf_bs_id: str = file_meta_dict[cbio_dx][cbio_tum_id]["kf_tum_id"]
            gene_fname: str = orig_fname + w_gene_suffix
            print(
                f"Found relevant cnv to process {file_meta_dict[cbio_dx][cbio_tum_id]['kf_tum_id']} {file_meta_dict[cbio_dx][cbio_tum_id]['kf_norm_id']} {gene_fname}",
                file=sys.stderr,
            )
            sys.stderr.flush()
            s_list.append(cbio_tum_id)
            if manifest is not None:
                ploidy: str = get_ploidy(info_dir + "/" + manifest.loc[kf_bs_id]["File_Name"])
                ploidy_dict[cbio_tum_id] = ploidy
            else:
                ploidy_dict[cbio_tum_id] = "2"
            cur_cnv_dict = process_cnv(cnv_dir + gene_fname, cur_cnv_dict, cbio_tum_id)
        new_cnv.write("Hugo_Symbol\t" + "\t".join(s_list) + "\n")
        for gene in cur_cnv_dict:
            new_cnv.write(gene)
            for samp in s_list:
                if samp in cur_cnv_dict[gene]:
                    new_cnv.write("\t" + cur_cnv_dict[gene][samp])
                else:
                    new_cnv.write("\t" + ploidy_dict[samp])
            new_cnv.write("\n")
        new_cnv.close()
        return 0, cbio_tum_id
    except Exception as e:
        print(e, file=sys.stderr)
        return 1, cbio_tum_id


parser = argparse.ArgumentParser(
    description="Merges files in gene <tab> entrez id <tab> copy number format into a genes-by-sample copy number table"
)
parser.add_argument(
    "-t",
    "--table",
    action="store",
    dest="table",
    help="Table with cbio project, kf bs ids, cbio IDs, and file names",
)
parser.add_argument(
    "-n",
    "--cnv-dir-gene",
    action="store",
    dest="cnv_dir",
    help="cnv as gene file directory",
)
parser.add_argument(
    "-i",
    "--info_dir",
    action="store",
    dest="info_dir",
    help="If info files available, provide here",
)
parser.add_argument(
    "-j",
    "--config",
    action="store",
    dest="config_file",
    help="json config file with data types and data locations",
)


args = parser.parse_args()
TOOL_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

with open(args.config_file) as f:
    config_data = json.load(f)
config_data: dict = resolve_config_paths(config_data, TOOL_DIR)

cnv_dir: str = args.cnv_dir
info_dir: str = args.info_dir
if cnv_dir[-1] != "/":
    cnv_dir += "/"
manifest = None
if info_dir is not None:
    manifest: Dataframe = pd.read_csv(args.table, sep="\t")
    manifest.set_index(["T_CL_BS_ID"], inplace=True)
    manifest = manifest.loc[manifest["File_Type"] == "info"]

else:
    sys.stderr.write("No info file given, will assume ploidy 2\n")
w_gene_suffix: str = ".CNVs.Genes.copy_number"
out_dir: str = "merged_cnvs/"
try:
    os.mkdir(out_dir)
except:
    sys.stderr.write("output dir already exists\n")
file_meta_dict: dict[str, dict[str, dict[str, str]]] = get_file_metadata(args.table, "cnv")
with concurrent.futures.ProcessPoolExecutor(config_data["cpus"]) as executor:
    results = {
        executor.submit(process_table, cbio_dx, file_meta_dict): cbio_dx
        for cbio_dx in file_meta_dict
    }
    for result in concurrent.futures.as_completed(results):
        if result.result()[0]:
            print(f"Failed processing {result.result()[1]}", file=sys.stderr)
            sys.exit(1)
