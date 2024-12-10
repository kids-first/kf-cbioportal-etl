""" 
Wrapper script to run all scripts in ETL pipeline

Run on Mgmt-Console-Dev-chopd3bprod@684194535433 EC2 instance

Run a specific step, or 2 steps, or all steps:
python3 run_pipeline.py --step [1, 2 3, or all]

Steps list:
1 - Get study metadata
2 - Get files from manifest
3 - Check downloaded files
4 - Build genomic file package

Example: 
pip install -e ~/tools/kf-cbioportal-etl/
cbio_etl_runner --help

cbio_etl_runner \
    --steps all \
    --db-ini /path/to/db.ini \
    --study oligo_nation \
    --meta-config-file oligo_nation_case_meta_config.json \
    --data-processing-config oligo_nation_data_processing_config.json \
    --sbg-profile default \
    --study-flag kf
"""

import argparse
import os
import subprocess
import sys

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(TOOL_DIR, "scripts")
CONFIGS_DIR = os.path.join(TOOL_DIR, "STUDY_CONFIGS")
REFS_DIR = os.path.join(TOOL_DIR, "REFS")

def resolve_config_path(provided_path, default_path):
    """
    Resolve the configuration path based on user input:
    1. If no value is provided, use the default (CONFIGS_DIR/{study}_{config_name}.json)
    2. If a file name is provided, check if it's in cwd or CONFIGS_DIR
    3. If an absolute path is provided, use it as is.
    """
    if not provided_path:
        return default_path
    elif os.path.isabs(provided_path):
        return provided_path
    else:
        cwd_path = os.path.join(os.getcwd(), provided_path)
        config_dir_path = os.path.join(CONFIGS_DIR, provided_path)
        if os.path.exists(cwd_path):
            return cwd_path
        elif os.path.exists(config_dir_path):
            return config_dir_path
        else:
            raise FileNotFoundError(
                f"Configuration file '{provided_path}' not found in CWD or CONFIGS_DIR."
            )

def setup_parser():
    parser = argparse.ArgumentParser(description="Run cBioPortal ETL pipeline")

    # Arguments for all steps
    parser.add_argument(
        "--steps",
        nargs="+",
        type=str,
        help="Steps to execute (e.g., 1 2 3 or all)",
        choices=[str(i) for i in range(1, 6)] + ["all"]
    )

    # Arguments for Step 1 - get_study_metadata.py
    parser.add_argument("-db", "--db-ini", required=True, help="Database config file")
    parser.add_argument("-p", "--profile", default="postgresql", help="Profile name (default: postgresql)")
    parser.add_argument("-mc", "--meta-config-file", required=False, help="Metadata configuration file. Default: value inputted for --study + '_case_meta_config.json'")
    parser.add_argument("-r", "--ref-dir", required=False, help="Reference directory. Defaults to tool's ref dir if not provided.")
    parser.add_argument("-a", "--all", required=False, action="store_true", help="Include all relevant files, not just status=active files, NOT RECOMMENDED")

    # Arguments for Step 2 - diff_studies.py
    parser.add_argument("-u", "--url", default="https://pedcbioportal.kidsfirstdrc.org/api/v2/api-docs", help="URL to search against")
    parser.add_argument("-s", "--study", required=True, help="Cancer study ID")
    parser.add_argument("-at", "--token", required=False, help="Token file obtained from Web API. Required if running Step 2")
    parser.add_argument("-ds", "--datasheet-sample", default=os.path.join(os.getcwd(), "datasheets/data_clinical_sample.txt"), help="File containing cBio-formatted sample metadata (default: datasheets/data_clinical_sample.txt from Step 1 output)")
    parser.add_argument("-dp", "--datasheet-patient", default=os.path.join(os.getcwd(), "datasheets/data_clinical_patient.txt"), help="File containing cBio-formatted patient metadata (default: datasheets/data_clinical_patient.txt from Step 1 output)")

    # Arguments for Step 3 - get_files_from_manifest.py
    parser.add_argument("-m", "--manifest", default="cbio_file_name_id.txt", help="Manifest file (default: cbio_file_name_id.txt from Step 1 output)")
    parser.add_argument("-f", "--file-types", required=False, help="Comma-separated file types to download")
    parser.add_argument("-t", "--aws-tbl", required=False, help="AWS table with bucket name and keys")
    parser.add_argument("-sp", "--sbg-profile", required=False, help="SBG profile name")
    parser.add_argument("-c", "--cbio-manifest", required=False, help="cBio manifest to limit downloads")
    parser.add_argument("-ao", "--active-only", required=False, action="store_true", help="Only include active files")
    parser.add_argument("-rm", "--rm-na", required=False, action="store_true", help="Remove entries where file_id and s3_path are NA")
    parser.add_argument("-d", "--debug", required=False, action="store_true", help="Enable debug mode")
    parser.add_argument("-o", "--overwrite", required=False, action="store_true", help="Overwrite files if they already exist")

    # Arguments for Step 5 - genomics_file_cbio_package_build.py
    parser.add_argument("-dpc", "--data-processing-config", required=False, help="Data processing configuration file. Default: value inputted for --study + '_data_processing_config.json'")
    parser.add_argument("-sf", "--study-flag", required=True, choices=["both", "kf", "dgd"], help="Study flag")
    parser.add_argument("-l", "--legacy", required=False, action="store_true", help="Enable legacy mode")

    return parser

def execute_step(step_number, args, steps):
    step = steps.get(step_number)
    if not step:
        sys.stderr.write(f"Step {step_number} not defined in the pipeline.\n")
        return

    sys.stderr.write(f"\nRunning Step {step_number}: {step['description']}\n\n")

    command = step["command"](args)
    try:
        subprocess.run(command, check=True)
        sys.stderr.write(f"\nCompleted Step {step_number}\n")
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"\nError in Step {step_number}. Ensure all dependencies are satisfied and try again.\n")
        sys.exit(1)
            
def main():
    parser = setup_parser()
    args = parser.parse_args()

    if args.ref_dir is None:
        args.ref_dir = REFS_DIR

    default_meta_config = os.path.join(CONFIGS_DIR, f"{args.study}_case_meta_config.json")
    default_data_config = os.path.join(CONFIGS_DIR, f"{args.study}_data_processing_config.json")
    args.meta_config_file = resolve_config_path(args.meta_config_file, default_meta_config)
    args.data_processing_config = resolve_config_path(args.data_processing_config, default_data_config)

    steps = {
        1: {
            "description": "Get study metadata",
            "output_files": [
                "datasheets/",
                "cbio_file_name_id.txt"
            ],
            "command": lambda args: (
                [
                    "python3",
                    os.path.join(SCRIPTS_DIR, "get_study_metadata.py"),
                    "-d", args.db_ini,
                    "-p", args.profile,
                    "-c", args.meta_config_file,
                    "-r", args.ref_dir
                ]
                + (["-a"] if args.all else [])
            )
        },
        2: {
            "description": "Compare studies on the server and local",
            "output_files": [
                "patient_portal_v_build.txt", 
                "sample_portal_v_build,txt"
            ], 
            "command": lambda args: (
                [
                    "python3",
                    os.path.join(SCRIPTS_DIR, "diff_studies.py"),
                    "-u", args.url,
                    "-c", args.study,
                    "-t", args.token,
                    "-s", args.datasheet_sample, 
                    "-p", args.datasheet_patient
                ]
            )
        },
        3: {
            "description": "Get files from manifest",
            "output_files": [
                "RSEM_gene/",
                "annofuse_filtered_fusions_tsv/",
                "annotated_public_outputs/",
                "ctrlfreec_bam_seg/",
                "ctrlfreec_info/",
                "manifest_subset.tsv"
            ],
            "command": lambda args: (
                [
                    "python3",
                    os.path.join(SCRIPTS_DIR, "get_files_from_manifest.py"),
                    "-m", args.manifest
                ]
                + (["-f", args.file_types] if args.file_types else [])
                + (["-t", args.aws_tbl] if args.aws_tbl else [])
                + (["-s", args.sbg_profile] if args.sbg_profile else [])
                + (["-c", args.cbio_manifest] if args.cbio_manifest else [])
                + (["-a"] if args.active_only else [])
                + (["-r"] if args.rm_na else [])
                + (["-d"] if args.debug else [])
                + (["-o"] if args.overwrite else [])
            )
        },
        4: {
        "description": "Check all intended files were downloaded",
        "output_files": ["missing_files.txt"],  
        "command": lambda args: (
            [
                "python3",
                os.path.join(SCRIPTS_DIR, "check_downloads.py"),
                "-m", "manifest_subset.tsv"
            ]
        )
        },
        5: {
            "description": "Build genomic file package",
            "output_files": [
                "MAFS/",
                "converted_cnvs/",
                "ctrlfreec_pval/",
                "merged_cnvs/",
                "merged_fusion/",
                "merged_mafs/",
                "merged_rsem/",
                "processed/",
                "cnv_gene_annot.log",
                "cnv_merge_seg.log",
                "collate_mafs.log",
                "convert_fusion_as_sv.log",
                "load_package_create.log",
                "merge_cnv_gene.log",
                "merge_cnv_gistic.log",
                "rna_merge_rename_expression.log",
                "validator.errs",
                "validator.out"
            ],
            "command": lambda args: (
                [
                    "python3",
                    os.path.join(SCRIPTS_DIR, "genomics_file_cbio_package_build.py"),
                    "-t", args.manifest,
                    "-c", args.meta_config_file,
                    "-d", args.data_processing_config,
                    "-f", args.study_flag
                ]
                + (["-l"] if args.legacy else [])
            )
        }
    }

    if "all" in args.steps:
        args.steps = list(steps.keys()) 
    args.steps = [int(step) for step in args.steps]

    for step_number in args.steps:
        execute_step(step_number, args, steps)


if __name__ == "__main__":
    main()
