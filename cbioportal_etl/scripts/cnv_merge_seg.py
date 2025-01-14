#!/usr/bin/env python3
import sys
import os
import argparse
import json
from get_file_metadata_helper import get_file_metadata
import concurrent.futures


def resolve_config_paths(config, tool_dir):
    """
    Resolve paths dynamically based on assumptions:
    - Paths starting with 'scripts/' or 'REFS/' are relative to the tool directory.
    """
    for key, value in config.items():
        if isinstance(value, dict):
            resolve_config_paths(value, tool_dir)
        elif isinstance(value, str) and value.startswith(("REFS/", "scripts/", "external_scripts/")):
            config[key] = os.path.abspath(os.path.join(tool_dir, value))

    return config


def process_seg(cur_seg, new_seg, cbio_tum_id):
    cur_seg = open(cur_seg)
    next(cur_seg)
    for seg in cur_seg:
        data = seg.rstrip("\n").split("\t")
        data[0] = cbio_tum_id
        # cbio does not like leading chr in chromosome name
        data[1] = data[1].replace("chr", "")
        new_seg.write("\t".join(data) + "\n")


def process_tbl(cbio_dx, file_meta_dict):
    try:
        x = 0
        # project/disease name should be name of directory hosting datasheet
        sys.stderr.write("Processing " + cbio_dx + " project" + "\n")
        new_seg = open(out_dir + cbio_dx + ".merged_seg.txt", "w")
        new_seg.write(seg_file_header)
        for cbio_tum_id in file_meta_dict[cbio_dx]:
            cbio_norm_id = file_meta_dict[cbio_dx][cbio_tum_id]["cbio_norm_id"]
            fname = file_meta_dict[cbio_dx][cbio_tum_id]["fname"]
            sys.stderr.write(
                "Found relevant seg to process for {} {} {} {} {}\n".format(
                    cbio_tum_id,
                    cbio_norm_id,
                    file_meta_dict[cbio_dx][cbio_tum_id]["kf_tum_id"],
                    file_meta_dict[cbio_dx][cbio_tum_id]["kf_norm_id"],
                    fname,
                )
            )
            sys.stderr.flush()
            process_seg(seg_dir + fname, new_seg, cbio_tum_id)
            x += 1
        sys.stderr.write(
            "Completed processing " + str(x) + " entries in " + cbio_dx + "\n"
        )
        new_seg.close()
    except Exception as e:
        sys.stderr.write(e)
        sys.exit()


parser = argparse.ArgumentParser(description="Merge seg files")
parser.add_argument(
    "-t",
    "--table",
    action="store",
    dest="table",
    help="Table with cbio project, kf bs ids, cbio IDs, and file names",
)
parser.add_argument(
    "-m", "--seg-dir", action="store", dest="seg_dir", help="seg file directory"
)
parser.add_argument(
    "-j",
    "--config",
    action="store",
    dest="config_file",
    help="json config file with data types and data locations",
)

args = parser.parse_args()
if __name__ == "__main__":
    TOOL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    with open(args.config_file) as f:
        config_data = json.load(f)
    config_data = resolve_config_paths(config_data, TOOL_DIR)
    
    seg_file_header = "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n"
    # get maf file ext
    seg_dir = args.seg_dir
    if seg_dir[-1] != "/":
        seg_dir += "/"
    file_meta_dict = get_file_metadata(args.table, "seg")

    out_dir = "merged_cnvs/"
    try:
        os.mkdir(out_dir)
    except:
        sys.stderr.write("output dir already exists\n")

    with concurrent.futures.ProcessPoolExecutor(config_data["cpus"]) as executor:
        results = {
            executor.submit(process_tbl, cbio_dx, file_meta_dict): cbio_dx
            for cbio_dx in file_meta_dict
        }

    sys.stderr.write("Done, check logs\n")
