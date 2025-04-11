#!/usr/bin/env python3
"""Genomics file + clinical data file to cBio package conversion workflow script.

It is designed to work with standard KF somatic workflow outputs as well as DGD outputs.
Clinical files should have been produced ahead of time, while supporting sample ID file manifest, case_meta config json file, and data json config.
See README for prerequisite details.
"""

import argparse
import json
import os
import subprocess
import sys
import time
from functools import partial

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths


def log_cmd(cmd):
    """Print output commands run to stderr."""
    print(cmd, file=sys.stderr)
    sys.stderr.flush()


def check_status(status, data_type, run_status):
    if status:
        print(f"Something when wrong while processing the {data_type} shutting down other running procs", file=sys.stderr)
        for key in run_status:
            if run_status[key] is None:
                run_status[key].kill()
        sys.exit(1)
    else:
        sys.stderr.write("Processing " + data_type + " successful!\n")
        sys.stderr.flush()


def process_maf(
    maf_loc_dict: dict[str, str | list[str]],
    cbio_id_table: str,
    data_config_file: str,
    dgd_status: str,
    script_dir: str,
    run_status: dict[str, subprocess.Popen],
) -> None:
    """Collate and process pbta/kf style mafs. Call append if also adding dgd data.

    Args:
        maf_loc_dict: Config dict entry with dbt file_type names and file header location
        cbio_id_table: Path of ETL table with IDs and file locations
        data_config_file: Path of file with ETL reference locations
        dgd_status: values of kf, both, and dgd to determine whether to append data or not
        script_dir: path of dir containing ETL processing scripts
        run_status: Dict to track subprocess call statuses

    """
    print("Processing maf files", file=sys.stderr)

    maf_header = maf_loc_dict["header"]
    for maf_type, maf_dir in maf_loc_dict.items():
        if maf_type in ["header", "dgd"]:
            continue
        if not maf_dir:
            print(f"Skipping {maf_type} as it is not defined in the config", file=sys.stderr)
            continue
        if isinstance(maf_dir, list):
            maf_dir = ",".join(maf_dir)
        maf_cmd = f"python3 {os.path.join(script_dir, 'maf_merge.py')} -t {cbio_id_table} -i {maf_header} -m {maf_dir} -j {data_config_file} -f {dgd_status} 2> collate_mafs.log"
        log_cmd(maf_cmd)
        run_status["maf"] = subprocess.Popen(maf_cmd, shell=True)


def process_append_dgd_maf(
    maf_loc_dict: dict[str, str | list[str]],
    cbio_id_table: str,
    cbio_study_id: str,
    script_dir: str,
) -> subprocess.Popen:
    """Append DGD mafs to existing collated kf maf.

    Args:
        maf_loc_dict: Config dict entry with dbt file_type names and file header location
        cbio_id_table: Path of ETL table with IDs and file locations
        cbio_study_id: Name of study being processed
        script_dir: path of dir containing ETL processing scripts

    Returns:
        Subprocess run of append python script

    """
    print("Appending DGD to exsting KF maf", file=sys.stderr)
    maf_header = maf_loc_dict["header"]
    in_maf_dir = maf_loc_dict["dgd"]
    append_maf = "merged_mafs/" + cbio_study_id + ".maf"
    append_cmd = f"python3 {os.path.join(script_dir, 'add_dgd_maf_to_pbta.py')} -i {maf_header} -m {in_maf_dir} -t {cbio_id_table} >> {append_maf} 2> dgd_append_maf.log"
    log_cmd(append_cmd)
    return subprocess.Popen(append_cmd, shell=True)


def process_cnv(
    cnv_loc_dict: dict[str, str],
    data_config_file: str,
    cbio_id_table: str,
    script_dir: str,
    run_status: dict[str, subprocess.Popen],
) -> None:
    """Add gene info to CNV calls, merge into table, and create GISTIC-style output.

    Args:
        cnv_loc_dict: Config dict entry with dbt file_type names and file header location
        cbio_id_table: Path of ETL table with IDs and file locations
        data_config_file: Path of file with ETL reference locations
        script_dir: path of dir containing ETL processing scripts
        run_status: Dict to track subprocess call statuses

    """
    sys.stderr.write("Processing CNV calls\n")
    cnv_dir: str = cnv_loc_dict["pval"]
    gene_annot_cmd = f"python3 {os.path.join(script_dir, 'cnv_1_genome2gene.py')} -d {cnv_dir} -j {data_config_file} 2> cnv_gene_annot.log"
    log_cmd(gene_annot_cmd)
    run_status["convert_to_gene"] = subprocess.Popen(gene_annot_cmd, shell=True)
    run_status["convert_to_gene"].wait()
    info_dir: str = cnv_loc_dict["info"]
    merge_gene_cmd = f"python3 {os.path.join(script_dir, 'cnv_2_merge.py')} -t {cbio_id_table} -n converted_cnvs -j {data_config_file}"
    if info_dir != "":
        merge_gene_cmd += f" -i {info_dir}"
    merge_gene_cmd += " 2> merge_cnv_gene.log"
    log_cmd(merge_gene_cmd)
    run_status["cnv_merge_gene"] = subprocess.Popen(merge_gene_cmd, shell=True)
    run_status["cnv_merge_gene"].wait()
    gistic_style_cmd = f"python3 {os.path.join(script_dir, 'cnv_3_gistic_style.py')} -d merged_cnvs -j {data_config_file} -t {cbio_id_table}"
    if info_dir != "":
        gistic_style_cmd += f" -i {info_dir}"
    gistic_style_cmd += " 2> merge_cnv_gistic.log"
    log_cmd(gistic_style_cmd)
    run_status["gistic"] = subprocess.Popen(gistic_style_cmd, shell=True)

    if "seg" in cnv_loc_dict:
        run_status["seg"] = process_seg(cnv_loc_dict, cbio_id_table, data_config_file, script_dir)


def process_seg(
    cnv_loc_dict: dict[str, str], cbio_id_table: str, data_config_file: str, script_dir: str
) -> subprocess.Popen:
    """Collate and process CNV seg files.

    Args:
        cnv_loc_dict: Config dict entry with dbt file_type names and file header location
        cbio_id_table: Path of ETL table with IDs and file locations
        data_config_file: Path of file with ETL reference locations
        script_dir: path of dir containing ETL processing scripts

    Returns:
        subprocess of process seg

    """
    sys.stderr.write("Processing CNV seg calls\n")
    seg_dir = cnv_loc_dict["seg"]
    merge_seg_cmd = f"python3 {os.path.join(script_dir, 'cnv_merge_seg.py')} -t {cbio_id_table} -m {seg_dir} -j {data_config_file} 2> cnv_merge_seg.log"
    log_cmd(merge_seg_cmd)
    return subprocess.Popen(merge_seg_cmd, shell=True)


def process_rsem(
    rsem_dir: str, cbio_id_table: str, script_dir: str, run_status: dict[str, subprocess.Popen]
) -> None:
    """Merge rsem results by FPKM, calculate z-scores.

    Args:
        rsem_dir: Path to rsem data
        cbio_id_table: Path of ETL table with IDs and file locations
        script_dir: path of dir containing ETL processing scripts
        run_status: Dict to track subprocess call statuses

    """
    print("Processing RNA expression data", file=sys.stderr)
    merge_rsem_cmd = f"python3 {os.path.join(script_dir, 'rna_merge_rename_expression.py')} -t {cbio_id_table} -r {rsem_dir} 2> rna_merge_rename_expression.log"
    log_cmd(merge_rsem_cmd)
    run_status["rsem_merge"] = subprocess.Popen(merge_rsem_cmd, shell=True)


def process_kf_fusion(
    fusion_dir: str,
    cbio_id_table: str,
    mode: str,
    script_dir: str,
    run_status: dict[str, subprocess.Popen],
) -> None:
    """Collate and process annoFuse output.

    Args:
        fusion_dir: Path to fusion data
        cbio_id_table: Path of ETL table with IDs and file locations
        mode: describe source - openX or kfprod or dgd
        script_dir: path of dir containing ETL processing scripts
        run_status: Dict to track subprocess call statuses

    """
    print("Processing KF fusion calls", file=sys.stderr)
    fusion_cmd = f"python3 {os.path.join(script_dir, 'convert_fusion_as_sv.py')} -t {cbio_id_table} -f {fusion_dir} -m {mode} 2> convert_fusion_as_sv.log"
    log_cmd(fusion_cmd)
    run_status["fusion"] = subprocess.Popen(fusion_cmd, shell=True)


def process_dgd_fusion(
    cbio_id_table: str, fusion_dir: str, dgd_status: str, script_dir: str, cbio_study_id: str
) -> subprocess.Popen:
    """Collate process DGD fusion output.

    Append if part of a KF/PBTA load, make solo file if not
    Args:
        fusion_dir: Path to fusion data
        cbio_id_table: Path of ETL table with IDs and file locations
        dgd_status: kf, both, dgd
        script_dir: path of dir containing ETL processing scripts

    Returns:
        Subprocess run of DGD fusion processing step

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
    return subprocess.Popen(dgd_fusion_cmd, shell=True)


def run_py(args):
    TOOL_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    with open(args.study_config) as f:
        config_data: dict = json.load(f)
    config_data = resolve_config_paths(config_data, TOOL_DIR)

    if args.dgd_status != "kf":
        try:
            fusion_file = config_data["file_loc_defs"]["dgd_fusion"]
            fusion_url = config_data["file_loc_defs"]["dgd_fusion_url"]
        except KeyError as e:
            raise KeyError(f"Missing required config key: {e}. Please ensure 'dgd_fusion' and 'dgd_fusion_url' are defined in file_loc_defs.")

        if not os.path.exists(fusion_file):
            sys.stderr.write(f"{fusion_file} not found. Downloading from {fusion_url}...\n")
            try:
                subprocess.run(
                    f"curl -sSL -o {fusion_file} {fusion_url}",
                    shell=True,
                    check=True
                )
                sys.stderr.write(f"Downloaded {fusion_file} successfully.\n")
            except subprocess.CalledProcessError as e:
                sys.stderr.write(f"Failed to download {fusion_file}: {e}\n")
                raise

    cbio_study_id: str = config_data["study"]["cancer_study_identifier"]
    # iterate through config file - file should only have keys related to data to be loaded
    script_dir: str = os.path.join(TOOL_DIR, config_data["script_dir"])
    # Store a list of run priorities to ensure historically slower jobs kick off first using partial:
    # https://stackoverflow.com/questions/59221490/can-you-store-functions-with-parameters-in-a-list-and-call-them-later-in-python
    run_priority: list[str] = ["rsem", "mafs", "fusion", "cnvs"]
    run_queue: dict[str, partial] = {}
    run_status: dict[str, subprocess.Popen] = {}
    for key in config_data:
        if key.startswith("merged_"):
            data_type: str = "_".join(key.split("_")[1:])
            if data_type == "mafs":
                run_queue["mafs"] = partial(
                    process_maf,
                    config_data["file_loc_defs"]["mafs"],
                    args.manifest,
                    args.study_config,
                    args.dgd_status,
                    script_dir,
                    run_status,
                )
            elif data_type == "rsem":
                run_queue["rsem"] = partial(
                    process_rsem,
                    config_data["file_loc_defs"]["rsem"],
                    args.manifest,
                    script_dir,
                    run_status,
                )
            elif data_type == "fusion":
                # Status both works for...both, only when one is specifically picked should one not be run
                if args.dgd_status != "dgd":
                    run_queue["fusion"] = partial(
                        process_kf_fusion,
                        config_data["file_loc_defs"]["fusion"],
                        args.manifest,
                        "kfprod",
                        script_dir,
                        run_status,
                    )

            elif data_type == "cnvs":
                run_queue["cnvs"] = partial(
                    process_cnv,
                    config_data["file_loc_defs"]["cnvs"],
                    args.study_config,
                    args.manifest,
                    script_dir,
                    run_status,
                )

    for job in run_priority:
        if job in run_queue:
            run_queue[job]()

    # Wait for concurrent processes to finish and report statuses
    x: int = 30
    n: int = 0
    done: bool = False

    # cheat flag to add append dgd maf once KF/PBTA maf is done when set to both
    dgd_fusion_append: int = 0
    while not done:
        time.sleep(x)
        n += x
        # don't need to check the same process over and over again if done
        rm_keys: list[str] = []
        done = True
        print(f"{n} seconds have passed", file=sys.stderr)
        sys.stderr.flush()

        if dgd_fusion_append:
            run_queue["dgd_fusion_append"] = partial(
                process_dgd_fusion,
                args.manifest,
                config_data["file_loc_defs"]["dgd_fusion"],
                args.dgd_status,
                script_dir,
                cbio_study_id,
            )
            run_queue["dgd_fusion_append"]()
            sys.stderr.write("dgd status was both, appending dgd fusions\n")
            # don't want to restart this job, let regular logic take over
            dgd_fusion_append = 0

        for data_type, process in run_status.items():
            process.poll()
            print(f"Checking {data_type} status", file=sys.stderr)
            if process.returncode is not None:
                check_status(process.returncode, data_type, run_status)
                rm_keys.append(data_type)
                if data_type == "fusion" and args.dgd_status == "both":
                    dgd_fusion_append = 1
                    done = False
            else:
                done = False
        for data_type in rm_keys:
            print(f"Removed {data_type} from check queue as task is complete", file=sys.stderr)
            del run_status[data_type]

    # Run final package builder script
    print("Creating load packages", file=sys.stderr)
    pck_cmd = f"python3 {os.path.join(script_dir, 'organize_upload_packages.py')} -o processed -c {args.study_config}"
    if args.add_data:
        pck_cmd += " -ad"
    pck_cmd += " 2> load_package_create.log"
    exit_status = subprocess.call(pck_cmd, shell=True)
    check_status(exit_status, "load package", run_status)

    # Run cbioportal data validator
    if not args.add_data:
        print("Validating load packages", file=sys.stderr)
        sys.stderr.flush()
        validate = f"{config_data['cbioportal_validator']}  -s processed/{cbio_study_id} -n -v 2> validator.errs > validator.out"
        exit_status = subprocess.call(validate, shell=True)
        if exit_status:
            print(
                f"Validator quit with status {exit_status}. Check validator.errs and validator.out for more info",
                file=sys.stderr,
            )


def main():
    parser = argparse.ArgumentParser(
        description="Download files (if needed), collate genomic files, organize load package."
    )
    parser.add_argument(
        "-m",
        "--manifest",
        action="store",
        dest="manifest",
        help="Download file manifest with cbio project, kf bs ids, cbio IDs, and file names",
    )
    parser.add_argument(
        "-sc", "--study-config", action="store", dest="study_config", help="cbio study config file"
    )
    parser.add_argument(
        "-dgd",
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
        "-ad",
        "--add-data",
        action="store_true",
        dest="add_data",
        help="Flag to skip validation when running for add_data directory",
    )

    args = parser.parse_args()
    run_py(args)


if __name__ == "__main__":
    main()
