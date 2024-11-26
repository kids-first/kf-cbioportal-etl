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
    --db-ini ../db.ini \
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

def setup_parser():
    parser = argparse.ArgumentParser(description="Run cBioPortal ETL pipeline")

    # Arguments for all steps
    parser.add_argument(
        "--steps",
        nargs="+",
        type=str,
        help="Steps to execute (e.g., 1 2 3 or all)",
        choices=[str(i) for i in range(1, 4)] + ["all"]
    )

    # Arguments for Step 1 - get_study_metadata.py
    parser.add_argument("-db", "--db-ini", required=True, help="Database config file")
    parser.add_argument("-p", "--profile", default="postgresql", help="Profile name (default: postgresql)")
    parser.add_argument("-mc", "--meta-config-file", required=True, help="Metadata configuration file. Defaults to tool's STUDY_CONFIGS dir")
    parser.add_argument("-r", "--ref-dir", default=None, required=False, help="Reference directory. Defaults to tool's ref dir if not provided.")
    parser.add_argument("-a", "--all", required=False, action="store_true", help="Include all relevant files, not just status=active files, NOT RECOMMENDED")

    # Arguments for Step 2 - get_files_from_manifest.py
    parser.add_argument("-m", "--manifest", default="cbio_file_name_id.txt", help="Manifest file (default: cbio_file_name_id.txt)")
    parser.add_argument("-f", "--file-types", required=False, help="Comma-separated file types to download")
    parser.add_argument("-t", "--aws-tbl", required=False, help="AWS table with bucket name and keys")
    parser.add_argument("-s", "--sbg-profile", required=False, help="SBG profile name")
    parser.add_argument("-c", "--cbio-manifest", required=False, help="cBio manifest to limit downloads")
    parser.add_argument("-ao", "--active-only", required=False, action="store_true", help="Only include active files")
    parser.add_argument("-rm", "--rm-na", required=False, action="store_true", help="Remove entries where file_id and s3_path are NA")
    parser.add_argument("-d", "--debug", required=False, action="store_true", help="Enable debug mode")
    parser.add_argument("-o", "--overwrite", required=False, action="store_true", help="Overwrite files if they already exist")

    # Arguments for Step 4 - genomics_file_cbio_package_build.py
    parser.add_argument("-dpc", "--data-processing-config", required=True, help="Data processing configuration file. Defaults to tool's STUDY_CONFIGS dir")
    parser.add_argument("-sf", "--study-flag", required=True, choices=["both", "kf", "dgd"], help="Study flag")
    parser.add_argument("-l", "--legacy", required=False, action="store_true", help="Enable legacy mode")

    return parser


def execute_step(step_number, args, steps):
    step = steps.get(step_number)
    if not step:
        sys.stderr.write(f"Step {step_number} not defined in the pipeline.\n")
        return

    command = step["command"](args)
    output_files = step["output_files"]

    # Check if output files from this step already exist
    missing_outputs = [f for f in output_files if not os.path.exists(f)]
    if missing_outputs:
        sys.stderr.write(f"\nRunning Step {step_number}: {step['description']}\n\n")
        # sys.stderr.write(f"Command: {' '.join(command)}\n\n")  # For debugging
        try:
            subprocess.run(command, check=True)
            sys.stderr.write(f"\nCompleted Step {step_number}\n")
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"\nError in Step {step_number}: {e}\n")
            sys.exit(1)
    # else:
    #     sys.stderr.write(f"\nSkipping Step {step_number}: Outputs already exist.\n")  # Not sure if we want to always rerun with fresh results, or to skip steps that have already been ran


def main():
    parser = setup_parser()
    args = parser.parse_args()

    args.meta_config_file = os.path.join(CONFIGS_DIR, args.meta_config_file)
    args.data_processing_config = os.path.join(CONFIGS_DIR, args.data_processing_config)
    if args.ref_dir is None:
        args.ref_dir = REFS_DIR

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
        3: {
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
        4: {
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
