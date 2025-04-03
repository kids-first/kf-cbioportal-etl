#!/usr/bin/env python3
"""Merge seg files and reformat for cBio."""

import argparse
import concurrent.futures
import json
import os
import sys
from typing import IO

from get_file_metadata_helper import get_file_metadata

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths


def process_seg(cur_seg_fname: str, new_seg: IO, cbio_tum_id: str) -> None:
    """Replace BS ID with cBio id and remove leading chr.

    Args:
        cur_seg_fname: Input seg file name
        new_seg: output file handle
        cbio_tume_id: cBio sample ID to use
    """
    with open(cur_seg_fname) as cur_seg:
        next(cur_seg)
        for seg in cur_seg:
            data: list = seg.rstrip("\n").split("\t")
            data[0] = cbio_tum_id
            # cbio does not like leading chr in chromosome name
            data[1] = data[1].replace("chr", "")
            print("\t".join(data), file=new_seg)


def process_tbl(cbio_dx: str, file_meta_dict: dict[str, dict[str, dict[str, str]]]) -> None:
    """Iterate through seg files in manifest.

    Create output file based on project
    Args:
        cbio_dx: cBio disease/project name
        file_meta_dict: Dict that has been subset by file type from ETL file
    """
    try:
        x: int = 0
        # project/disease name should be name of directory hosting datasheet
        print(f"Processing {cbio_dx} project", file=sys.stderr)
        with open(f"{out_dir}{cbio_dx}.merged_seg.txt", "w") as new_seg:
            print(seg_file_header, file=new_seg)
            for cbio_tum_id in file_meta_dict[cbio_dx]:
                cbio_norm_id: str = file_meta_dict[cbio_dx][cbio_tum_id]["cbio_norm_id"]
                fname: str = file_meta_dict[cbio_dx][cbio_tum_id]["fname"]
                print(
                    f"Found relevant seg to process for {cbio_tum_id} {cbio_norm_id} {file_meta_dict[cbio_dx][cbio_tum_id]['kf_tum_id']} {file_meta_dict[cbio_dx][cbio_tum_id]['kf_norm_id']} {fname}",
                    file=sys.stderr,
                )
                sys.stderr.flush()
                process_seg(seg_dir + fname, new_seg, cbio_tum_id)
                x += 1
            print(f"Completed processing {x} entries in {cbio_dx}", file=sys.stderr)
    except Exception as e:
        print(f"{e}", file=sys.stderr)
        sys.exit()


parser = argparse.ArgumentParser(description="Merge seg files")
parser.add_argument(
    "-t",
    "--table",
    action="store",
    dest="table",
    help="Table with cbio project, kf bs ids, cbio IDs, and file names",
)
parser.add_argument("-m", "--seg-dir", action="store", dest="seg_dir", help="seg file directory")
parser.add_argument(
    "-j",
    "--config",
    action="store",
    dest="config_file",
    help="json config file with data types and data locations",
)

args = parser.parse_args()
if __name__ == "__main__":
    TOOL_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    with open(args.config_file) as f:
        config_data: dict = json.load(f)
    config_data = resolve_config_paths(config_data, TOOL_DIR)

    seg_file_header: str = "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean"
    # get maf file ext
    seg_dir: str = args.seg_dir
    if seg_dir[-1] != "/":
        seg_dir += "/"
    file_meta_dict: dict[str, dict[str, dict[str, str]]] = get_file_metadata(args.table, "seg")

    out_dir: str = "merged_cnvs/"
    os.makedirs(out_dir, exist_ok=True)

    with concurrent.futures.ProcessPoolExecutor(config_data["cpus"]) as executor:
        results: dict[concurrent.futures.Future[None], str] = {
            executor.submit(process_tbl, cbio_dx, file_meta_dict): cbio_dx
            for cbio_dx in file_meta_dict
        }

    print("Done, check logs", file=sys.stderr)
