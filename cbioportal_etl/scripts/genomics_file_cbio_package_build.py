#!/usr/bin/env python3
"""
This is a genomics file + clinical data file to cBio package conversion workflow script.
It is designed to work with standard KF somatic workflow outputs as well as DGD outputs.
Clinical files should have been produced ahead of time, while supporting sample ID file manifest, case_meta config json file, and data json config.
If there are fusion files, they also need a table outlining sequencing center location.
See README for prequisite details.
"""
import os
import sys
import argparse
import json
import subprocess
import time
from functools import partial
import pdb

def resolve_config_paths(config, tool_dir, cwd):
    """
    Resolve paths dynamically based on assumptions:
    - Paths starting with 'scripts/' or 'REFS/' are relative to the tool directory.
    - Other paths are assumed to be relative to the current working directory.
    """
    for key, value in config.items():
        if isinstance(value, dict):
            resolve_config_paths(value, tool_dir, cwd)
        elif isinstance(value, str) and not value.startswith(("http://", "https://")):
            if value.startswith(("REFS/", "scripts/", "external_scripts/")):
                config[key] = os.path.abspath(os.path.join(tool_dir, value))
            else:
                config[key] = os.path.abspath(os.path.join(cwd, value))
    return config

def log_cmd(cmd):
    """
    A silly litte helper function to output commands run to stderr
    """
    sys.stderr.write(cmd + "\n")
    sys.stderr.flush()

def check_status(status, data_type):
    if status:
        sys.stderr.write(
            "Something when wrong while processing the "
            + data_type + " shutting down other running procs\n"
        )
        for key in run_status:
            if run_status[key] is None:
                run_status[key].kill()
        exit(1)
    else:
        sys.stderr.write('Processing ' + data_type + " successful!\n")
        sys.stderr.flush()


def process_maf(maf_loc_dict, cbio_id_table, data_config_file, dgd_status):
    """
    Collate and process pbta/kf style mafs. Call append if also adding dgd data
    """
    sys.stderr.write("Processing maf files\n")
    maf_dir = maf_loc_dict["kf"]
    if args.dgd_status == "dgd":
        maf_dir = maf_loc_dict["dgd"]
    else:
        # KF can be in multiple palces
        maf_dir = ",".join(maf_dir)
    maf_header = maf_loc_dict["header"]
    maf_cmd = f"python3 {os.path.join(script_dir, 'maf_merge.py')} -t {cbio_id_table} -i {maf_header} -m {maf_dir} -j {data_config_file} -f {dgd_status} 2> collate_mafs.log"
    log_cmd(maf_cmd)
    run_status["maf"] = subprocess.Popen(maf_cmd, shell=True)


def process_append_dgd_maf(maf_loc_dict, cbio_id_table):
    """
    Append DGD mafs to existing collated kf maf
    """
    sys.stderr.write("Appending DGD to exsting KF maf\n")
    maf_header = maf_loc_dict["header"]
    in_maf_dir = maf_loc_dict["dgd"]
    append_maf = "merged_mafs/" + cbio_study_id + ".maf"
    append_cmd = (
        f"python3 {os.path.join(script_dir, 'add_dgd_maf_to_pbta.py')} -i {maf_header} -m {in_maf_dir} -t {cbio_id_table} >> {append_maf} 2> dgd_append_maf.log"
    )
    log_cmd(append_cmd)
    status = 0
    status = subprocess.Popen(append_cmd, shell=True)
    return status


def process_cnv(cnv_loc_dict, data_config_file, cbio_id_table):
    """
    Add gene info to CNV calls, merge into table, and create GISTIC-style output
    """
    sys.stderr.write("Processing CNV calls\n")
    cnv_dir = cnv_loc_dict["pval"]
    gene_annot_cmd = f"python3 {os.path.join(script_dir, 'cnv_1_genome2gene.py')} -d {cnv_dir} -j {data_config_file} 2> cnv_gene_annot.log"
    log_cmd(gene_annot_cmd)
    run_status["convert_to_gene"] = subprocess.Popen(gene_annot_cmd, shell=True)
    run_status["convert_to_gene"].wait()
    info_dir = cnv_loc_dict["info"]
    merge_gene_cmd = f"python3 {os.path.join(script_dir, 'cnv_2_merge.py')} -t {cbio_id_table} -n converted_cnvs -j {data_config_file}"
    if info_dir != "":
        merge_gene_cmd += " -i " + info_dir
    merge_gene_cmd += " 2> merge_cnv_gene.log"
    log_cmd(merge_gene_cmd)
    run_status["cnv_merge_gene"] = subprocess.Popen(merge_gene_cmd, shell=True)
    run_status["cnv_merge_gene"].wait()
    gistic_style_cmd = f"python3 {os.path.join(script_dir, 'cnv_3_gistic_style.py')} -d merged_cnvs -j {data_config_file} -t {cbio_id_table}"
    if info_dir != "":
        gistic_style_cmd += " -i " + info_dir
    gistic_style_cmd += " 2> merge_cnv_gistic.log"
    log_cmd(gistic_style_cmd)
    run_status["gistic"] = subprocess.Popen(gistic_style_cmd, shell=True)

    if "seg" in cnv_loc_dict:
        run_status["seg"] = process_seg(cnv_loc_dict, cbio_id_table, data_config_file)


def process_seg(cnv_loc_dict, cbio_id_table, data_config_file):
    """
    Collate and process CNV seg files
    """
    sys.stderr.write("Processing CNV seg calls\n")
    seg_dir = cnv_loc_dict["seg"]
    merge_seg_cmd = f"python3 {os.path.join(script_dir, 'cnv_merge_seg.py')} -t {cbio_id_table} -m {seg_dir} -j {data_config_file} 2> cnv_merge_seg.log"
    log_cmd(merge_seg_cmd)
    status = subprocess.Popen(merge_seg_cmd, shell=True)
    return status


def process_rsem(rsem_dir, cbio_id_table):
    """
    Merge rsem results by FPKM, calculate z-scores
    """
    sys.stderr.write("Processing RNA expression data\n")
    merge_rsem_cmd = f"python3 {os.path.join(script_dir, 'rna_merge_rename_expression.py')} -t {cbio_id_table} -r {rsem_dir} 2> rna_merge_rename_expression.log"
    log_cmd(merge_rsem_cmd)
    run_status["rsem_merge"] = subprocess.Popen(merge_rsem_cmd, shell=True)


def process_kf_fusion(fusion_dir, cbio_id_table, mode):
    """
    Collate and process annoFuse output
    """
    sys.stderr.write("Processing KF fusion calls\n")
    fusion_cmd = f"python3 {os.path.join(script_dir, 'convert_fusion_as_sv.py')} -t {cbio_id_table} -f {fusion_dir} -m {mode} 2> convert_fusion_as_sv.log"
    log_cmd(fusion_cmd)
    run_status["fusion"] = subprocess.Popen(fusion_cmd, shell=True)


def process_kf_fusion_legacy(fusion_dir, cbio_id_table, sq_file):
    """
    Collate and process annoFuse output using deprecated format
    """
    sys.stderr.write("Processing KF fusion calls\n")
    fusion_cmd = f"python3 {os.path.join(script_dir, 'rna_convert_fusion.py')} -t {cbio_id_table} -f {fusion_dir} -m annofuse -s {sq_file} 2> rna_convert_fusion.log"
    log_cmd(fusion_cmd)
    run_status["fusion_legacy"] = subprocess.Popen(fusion_cmd, shell=True)


def process_dgd_fusion(cbio_id_table, fusion_dir, dgd_status):
    """
    Collate process DGD fusion output
    """
    dgd_fusion_cmd = f"python3 {os.path.join(script_dir, 'convert_fusion_as_sv.py')} -t {cbio_id_table} -f {fusion_dir} -m dgd"
    if dgd_status == "both":
        sys.stderr.write("Appending DGD fusion calls\n")
        append_fusion = "merged_fusion/" + cbio_study_id + ".fusions.txt"
        dgd_fusion_cmd += " -a >> " + append_fusion
    else:
        sys.stderr.write("Processing DGD fusion calls\n")
    dgd_fusion_cmd += " 2> add_dgd_fusion.log"
    log_cmd(dgd_fusion_cmd)
    status = subprocess.Popen(dgd_fusion_cmd, shell=True)
    return status


# This part is commented out for now as there is no good solution yet for files stored in different account - need to be able to run chopaws to get key
# if(args.manifest != None):
#     sys.stderr.write('Download manifest given. Downloading files\n')
#     dl_file_cmd = config_data['script_dir'] + 'get_files_from_manifest.py -m ' + args.manifest + ' -f ' + join(config_data['dl_file_type_list'] + ' -p saml')
#     subprocess.Popen(dl_file_cmd, shell=True)

parser = argparse.ArgumentParser(
    description="Download files (if needed), collate genomic files, organize load package."
)
parser.add_argument(
    "-t",
    "--table",
    action="store",
    dest="table",
    help="Table with cbio project, kf bs ids, cbio IDs, and file names",
)
parser.add_argument(
    "-m",
    "--manifest",
    action="store",
    dest="manifest",
    help="Download file manifest, if needed",
)
parser.add_argument(
    "-c",
    "--cbio-config",
    action="store",
    dest="cbio_config",
    help="cbio case and meta config file",
)
parser.add_argument(
    "-d",
    "--data-config",
    action="store",
    dest="data_config",
    help="json config file with data types and data locations",
)
parser.add_argument(
    "-f",
    "--dgd-status",
    action="store",
    dest="dgd_status",
    help="Flag to determine load will have pbta/kf + dgd(both), kf/pbta only(kf), dgd-only(dgd)",
    default="both",
    const="both",
    nargs="?",
    choices=["both", "kf", "dgd"],
)
parser.add_argument(
    "-l",
    "--legacy", 
    default=False,
    action="store_true",
    dest="legacy",
    help="If set, will run legacy fusion output",
)

TOOL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CWD = os.getcwd()

args = parser.parse_args()
with open(args.data_config) as f:
    config_data = json.load(f)
config_data = resolve_config_paths(config_data, TOOL_DIR, CWD)
with open(args.cbio_config) as f:
    config_meta_case = json.load(f)

cbio_study_id = config_meta_case["study"]["cancer_study_identifier"]
# iterate through config file - file should only have keys related to data to be loaded
script_dir = os.path.join(TOOL_DIR, config_data["script_dir"])
run_status = {}
# Store a list of run priorities to ensure historically slower jobs kick off first using partial:
# https://stackoverflow.com/questions/59221490/can-you-store-functions-with-parameters-in-a-list-and-call-them-later-in-python
run_priority = ["rsem", "mafs", "fusion", "cnvs"]
run_queue = {}
for key in config_meta_case:
    if key.startswith("merged_"):
        data_type = "_".join(key.split("_")[1:])
        if data_type == "mafs":
            run_queue["mafs"] = partial(process_maf,
                config_data["file_loc_defs"]["mafs"],
                args.table,
                args.data_config,
                args.dgd_status
            )
        elif data_type == "rsem":
            run_queue["rsem"] = partial(process_rsem, config_data["file_loc_defs"]["rsem"], args.table)
        elif data_type == "fusion":
            # Status both works for...both, only when one is specifically picked should one not be run
            if args.dgd_status != "dgd":
                if not args.legacy:
                    run_queue["fusion"]=partial(process_kf_fusion,
                        config_data["file_loc_defs"]["fusion"],
                        args.table,
                        "kfprod"
                        )
                else:
                    run_queue["fusion"] = partial(process_kf_fusion_legacy, 
                        config_data["file_loc_defs"]["fusion"],
                        args.table,
                        config_data["file_loc_defs"]["fusion_sq_file"],
                    )
        elif data_type == "cnvs":
            run_queue["cnvs"] = partial(process_cnv,
                config_data["file_loc_defs"]["cnvs"], args.data_config, args.table
            )
for job in run_priority:
    if job in run_queue:
        run_queue[job]()

# Wait for concurrent processes to finish and report statuses
x = 30
n = 0
done = False

# cheat flag to add append dgd maf once KF/PBTA maf is done when set to both
dgd_maf_append = 0
dgd_fusion_append = 0
while not done:
    time.sleep(x)
    n += x
    # don't need to check the same process over and over again if done
    rm_keys = []
    done = True
    sys.stderr.write(str(n) + " seconds have passed\n")
    sys.stderr.flush()
    if dgd_maf_append:
        run_queue["dgd_maf_append"] = partial(process_append_dgd_maf, config_data["file_loc_defs"]["mafs"], args.table)
        run_queue["dgd_maf_append"]()
        sys.stderr.write("dgd status was both, appending dgd maf\n")
        # don't want to restart this job, let regular logic take over
        dgd_maf_append = 0
    if dgd_fusion_append:
        run_queue["dgd_fusion_append"] = partial(process_dgd_fusion, args.table, config_data["file_loc_defs"]["dgd_fusion"], args.dgd_status)
        run_queue["dgd_fusion_append"]()
        sys.stderr.write("dgd status was both, appending dgd fusions\n")
        # don't want to restart this job, let regular logic take over
        dgd_fusion_append = 0

    for key in run_status:
        run_status[key].poll()
        sys.stderr.write("Checking " + key + " status\n")
        if run_status[key].returncode != None:
            check_status(run_status[key].returncode, key)
            rm_keys.append(key)
            # not all studies that have DGD fusion have maf
            if key == 'maf' and args.dgd_status=='both' and 'dgd' in config_data["file_loc_defs"]["mafs"]:
                dgd_maf_append = 1
                done = False
            elif key == 'fusion' and args.dgd_status=='both':
                dgd_fusion_append = 1
                done = False

        else:
            done = False
    for key in rm_keys:
        sys.stderr.write("Removed " + key + " from check queue as task is complete\n")
        del run_status[key]



# Run final package builder script
sys.stderr.write("Creating load packages\n")
pck_cmd = f"python3 {os.path.join(script_dir, 'organize_upload_packages.py')} -o processed -c {args.cbio_config} 2> load_package_create.log"
exit_status = subprocess.call(pck_cmd, shell=True)
check_status(exit_status, "load package")

# Run cbioportal data validator
sys.stderr.write("Validating load packages\n")
sys.stderr.flush()
validate = (
    config_data["cbioportal_validator"]
    + " -s processed/"
    + cbio_study_id
    + " -n -v 2> validator.errs > validator.out"
)
exit_status = subprocess.call(validate, shell=True)
if exit_status:
    sys.stderr.write(
        "Validator quit with status "
        + str(exit_status)
        + ". Check validator.errs and validator.out for more info\n"
    )