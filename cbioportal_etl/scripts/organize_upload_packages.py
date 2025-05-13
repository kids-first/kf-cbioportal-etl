#!/usr/bin/env python3
"""Organize data in final upload format.

This includes creating soft links to all data files as well as meta files to describe data file
and case lists to describe number of samples of each data type
"""

import argparse
import json
import os
import subprocess
import sys

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths


def process_meta_study(meta_data: dict, output_dir: str) -> str:
    """Create meta_ format files for study cBio import.

    Args:
        meta_data: Pointer to dict within config dict with genomic file info
        output_dir: Location to create each meta file

    """
    study: str = meta_data["cancer_study_identifier"]
    # fields where key matches and value is just value of dict entry
    matched_fields: list[str] = [
        "type_of_cancer",
        "short_name",
        "reference_genome",
        "description",
        "groups",
        "citation",
        "pmid",
    ]
    with open(f"{output_dir}meta_study.txt", "w") as meta_study:
        print(f"cancer_study_identifier: {study}", file=meta_study)
        print(f"name: {meta_data['display_name']}", file=meta_study)
        for key in matched_fields:
            if key in meta_data:
                print(f"{key}: {meta_data[key]}", file=meta_study)
    return study


def process_meta_data(meta_data: dict, output_dir: str, canc_study_id: str) -> None:
    """Create meta_ format files for genomic data cBio import.

    Args:
        meta_data: Pointer to dict within config dict with genomic file info
        output_dir: Location to create each meta file
        canc_study_id: Name of study/project

    """
    try:
        for dtype in meta_data["dtypes"]:
            # pointer for easier readability and dict key traciing
            cur_data: dict = meta_data["dtypes"][dtype]

            # Check which expression file (FPKM or TPM) exists
            if dtype.startswith("counts_"):
                data_path = os.path.join(meta_data["dir"], f"{canc_study_id}.{cur_data['ext']}")
                if not os.path.exists(data_path):
                    print(f"Skipping {dtype} - {data_path} not found", file=sys.stderr)
                    continue

            # Skip healthy z-score if file does not exist
            if dtype.startswith("zscore_"):
                data_path = os.path.join(meta_data["dir"], f"{canc_study_id}.{cur_data['ext']}")
                if not os.path.exists(data_path):
                    print(f"Skipping {dtype} - {data_path} not found", file=sys.stderr)
                    continue

            cbio_name: str = cur_data["cbio_name"]
            parts: list[str] = cbio_name.split("_")
            attr_dict: dict = cur_data["meta_file_attr"]
            meta_name: str = "meta_" + "_".join(parts[1:])
            if attr_dict["datatype"] == "FUSION":
                meta_name = "meta_" + attr_dict["datatype"] + ".txt"
            with open(f"{output_dir}{meta_name}", "w") as meta_data_file:
                print(f"cancer_study_identifier: {canc_study_id}", file=meta_data_file)
                for meta_key, attribtute in attr_dict.items():
                    print(f"{meta_key}: {attribtute}", file=meta_data_file)
                print(f"data_filename: {cbio_name}", file=meta_data_file)
            # create data_ links to data
            cmd = f"ln {cwd}{meta_data['dir']}/{canc_study_id}.{cur_data['ext']} {output_dir}{cbio_name}"
            subprocess.call(cmd, shell=True)
    except Exception as e:
        print(f"{e} failed processing meta data file", file=sys.stderr)
        sys.exit(1)


def process_clinical_data(
    meta_data: dict, output_dir: str, canc_study_id: str, add_data_mode: bool = False
) -> None:
    """Create meta_ format files for clinical data cBio import.

    Args:
        meta_data: Pointer to dict within config dict with clinical file info
        output_dir: Location to create each meta file
        canc_study_id: Name of study/project
        add_data_mode: Flag whether creating in add data mode versus whole study

    """
    datasheets_dir: str = meta_data["dir"]
    if add_data_mode:
        existing_files = set(os.listdir(datasheets_dir))

    for dtype in meta_data["dtypes"]:
        try:
            cur_data: dict = meta_data["dtypes"][dtype]
            cbio_name: str = cur_data["cbio_name"]

            if add_data_mode and cbio_name not in existing_files:
                print(
                    f"Skipping {cbio_name} as it does not exist in {datasheets_dir}",
                    file=sys.stderr,
                )
                continue
            parts: list[str] = cbio_name.split("_")
            meta_name: str = "meta_" + "_".join(parts[1:])
            with open(output_dir + meta_name, "w") as meta_data_file:
                print(f"cancer_study_identifier: {canc_study_id}", file=meta_data_file)
                attr_dict: dict = cur_data["meta_file_attr"]
                for meta_key, attribute in attr_dict.items():
                    print(f"{meta_key}: {attribute}", file=meta_data_file)
                print(f"data_filename: {cbio_name}", file=meta_data_file)

            cmd = f"ln {cwd}{meta_data['dir']}/{cbio_name} {output_dir}{cbio_name}"
            subprocess.call(cmd, shell=True)
        except Exception as e:
            print(str(e) + " failed processing meta data file", file=sys.stderr)
            sys.exit(1)


def write_case_list(
    case_key: str, attr_dict: dict[str, str], sample_list: list[str], case_dir: str
) -> None:
    """Write case lists based on data type being described.

    Args:
        case_key: Data type
        attr_dict: Attribute dict dor that data type
        sample_list: List of samples relevant to data type
        case_dir: Output dir location

    """
    try:
        with open(case_dir + case_key + ".txt", "w") as case_file:
            print(f"cancer_study_identifier: {canc_study_id}", file=case_file)

            for key, value in attr_dict.items():
                if key == "case_list_description":
                    to_write = f"{key}: {value} ({len(sample_list)})"
                elif key == "stable_id":
                    to_write = f"{key}: {canc_study_id}_{value}"
                else:
                    to_write = f"{key}: {value}"
                print(to_write, file=case_file)
            joined_samples = '\t'.join(sample_list)
            print(f"case_list_ids: {joined_samples}", file=case_file)
    except Exception as e:
        print(f"{e}\nError writing case list for {case_key}", file=sys.stderr)


def create_case_lists(data_dict: dict[str, int], output_dir: str):
    """Iterate through config file for case list creation.

    Determine data types available, and initialize relevant sample lists for each data type and data type combo
    Args:
        data_dict: Dict with flags indicating data type present
    """
    try:
        case_dir = f"{output_dir}case_lists/"
        os.makedirs(case_dir, exist_ok=True)

        muts_list: list[str] = []
        cna_list: list[str] = []
        rna_list: list[str] = []
        fusion_list: list[str] = []

        if data_dict["merged_mafs"]:
            # Get samples from merged maf file
            muts_fname: str = (
                f"{output_dir}{config_data['merged_mafs']['dtypes']['mutation']['cbio_name']}"
            )
            with open(muts_fname) as muts_file:
                head: str = next(muts_file)
                head = next(muts_file)
                header: list[str] = head.rstrip("\n").split("\t")
                s_idx: int = header.index("Tumor_Sample_Barcode")
                for line in muts_file:
                    data: list[str] = line.rstrip("\n").split("\t")
                    muts_list.append(data[s_idx])
            muts_list = [*{*muts_list}]
            write_case_list("cases_sequenced", config_data["cases_sequenced"], muts_list, case_dir)
        # Initialize all cases with muts list, even if empty
        all_cases: list[str] = muts_list
        if data_dict["merged_cnvs"]:
            # Get samples from merged cnv file
            cna_fname: str = f"{output_dir}{config_data['merged_cnvs']['dtypes']['linear']['cbio_name']}"
            with open(cna_fname) as cna_file:
                head: str = next(cna_file)
                # assumes header is Hugo_symbols\tsample_name1\tsamplename2 etc, if entrez ID, will need to change!
                cna_list: list[str] = head.rstrip("\n").split("\t")[1:]
            write_case_list("cases_cna", config_data["cases_cna"], cna_list, case_dir)
            all_cases += cna_list
            # Create case list for samples with muts and cnv data
            muts_plus_cna: list[str] = list(set(muts_list) & set(cna_list))
            write_case_list("cases_cnaseq", config_data["cases_cnaseq"], muts_plus_cna, case_dir)

        if data_dict["merged_rsem"]:
            # Get samples from merged rsem file - is a simple gene by sample table
            rsem_dtypes = config_data["merged_rsem"]["dtypes"]
            counts_key = next((k for k in rsem_dtypes if k.startswith("counts_")), None)
            if counts_key:
                rna_fname = f"{output_dir}{rsem_dtypes[counts_key]['cbio_name']}"
            else:
                print("No counts_fpkm or counts_tpm found in merged_rsem", file=sys.stderr)
                sys.exit(1)
            with open(rna_fname) as rna_file:
                head: str = next(rna_file)
                rna_list: list[str] = head.rstrip("\n").split("\t")[1:]
            write_case_list(
                "cases_RNA_Seq_v2_mRNA", config_data["cases_RNA_Seq_v2_mRNA"], rna_list, case_dir
            )
            # loading mutations is a minimum, so if cna exists...3 way file can be made
            if len(cna_list) > 0:
                three_way: list[str] = list(set(muts_list) & set(cna_list) & set(rna_list))
                write_case_list(
                    "cases_3way_complete", config_data["cases_3way_complete"], three_way, case_dir
                )
            all_cases += rna_list

        if "merged_fusion" in data_keys and data_keys["merged_fusion"]:
            # Get samples from merge fusion file
            fusion_fname: str = (
                f"{output_dir}{config_data['merged_fusion']['dtypes']['fusion']['cbio_name']}"
            )
            with open(fusion_fname) as fusion_file:
                head: str = next(fusion_file)
                header: list[str] = head.rstrip("\n").split("\t")
                sample_id_colname: str = "Sample_Id"
                s_idx: int = header.index(sample_id_colname)
                for line in fusion_file:
                    data: list[str] = line.rstrip("\n").split("\t")
                    fusion_list.append(data[s_idx])
            fusion_list = [*{*fusion_list}]
            write_case_list("cases_sv", config_data["cases_sv"], fusion_list, case_dir)
            all_cases += fusion_list

        all_cases = [*{*all_cases}]
        write_case_list("cases_all", config_data["cases_all"], all_cases, case_dir)
    except Exception as e:
        print(f"{e}\nError creating case lists", file=sys.stderr)
        sys.exit(1)


parser = argparse.ArgumentParser(
    description="Create cases lists, meta files, and organize data for cbio upload."
    " It is assumed you are at the dir level of all input data files"
)
parser.add_argument(
    "-o", "--output_dir", action="store", dest="out_dir", help="output directory name"
)
parser.add_argument(
    "-c",
    "--config",
    action="store",
    dest="config_file",
    help="json config file with meta information; see REFS/case_meta_config.json example",
)
parser.add_argument(
    "-ad",
    "--add-data",
    action="store_true",
    dest="add_data",
    help="Flag to skip validation when running for add_data directory",
)

args = parser.parse_args()

cwd: str = os.getcwd() + "/"
TOOL_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

with open(args.config_file) as f:
    config_data: dict = json.load(f)
config_data: dict = resolve_config_paths(config_data, TOOL_DIR)

out_dir: str = args.out_dir if args.out_dir[-1] == "/" else args.out_dir + "/"
os.makedirs(out_dir, exist_ok=True)

try:
    study_id: str = config_data["study"]["cancer_study_identifier"]
    if os.path.isdir(config_data["data_sheets"]["dir"]):
        cur_dir: str = f"{out_dir}{config_data['study']['cancer_study_identifier']}/"
        os.makedirs(cur_dir, exist_ok=True)
        if not args.add_data:
            print(f"Creating meta study file for {study_id}", file=sys.stderr)
            canc_study_id = process_meta_study(config_data["study"], cur_dir)
        else:
            canc_study_id = study_id
        # track whether data type in study
        data_keys: dict[str, int] = {
            "merged_mafs": 0,
            "merged_cnvs": 0,
            "merged_rsem": 0,
            "merged_fusion": 0,
        }
        for key in data_keys:
            if key in config_data and config_data[key]["dir"] != "":
                data_keys[key] = 1
                process_meta_data(config_data[key], cur_dir, canc_study_id)
                print(f"Creating meta data files and links for {key}", file=sys.stderr)
            else:
                print(f"Skipping meta files for {key}", file=sys.stderr)
        print("Creating clinical meta sheets and link", file=sys.stderr)
        process_clinical_data(config_data["data_sheets"], cur_dir, canc_study_id, args.add_data)
        if not args.add_data:
            create_case_lists(data_keys, cur_dir)
    else:
        print(f"No datasheets for {study_id}, skipping!", file=sys.stderr)
except Exception as e:
    print(f"{e}\nerror processing files for {study_id}!", file=sys.stderr)
