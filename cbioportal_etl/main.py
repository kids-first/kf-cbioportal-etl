""" 
Wrapper script to run all scripts in ETL pipeline

Run on Mgmt-Console-Dev-chopd3bprod@684194535433 EC2 instance

Steps list:
1 - Get study metadata
2 - Compare studies on the server and local
3 - Get files from manifest
4 - Check downloaded files
5 - Build genomic file package

Example Usage: 
pip install /path/to/kf-cbioportal-etl/
cbioportal_etl \
    --steps all \
    --db-ini /path/to/db.ini \
    --study oligo_nation \
    --sbg-profile default 
"""

import argparse
import os
import sys
import subprocess
from cbioportal_etl.scripts.get_study_metadata import run_py as get_study_metadata
from cbioportal_etl.scripts.diff_studies import run_py as diff_studies
from cbioportal_etl.scripts.get_files_from_manifest import run_py as get_files_from_manifest
from cbioportal_etl.scripts.check_downloads import run_py as check_downloads
from cbioportal_etl.scripts.genomics_file_cbio_package_build import run_py as genomics_file_cbio_package_build


def fetch_validator_scripts(tool_dir):
    repo_url = "https://github.com/cBioPortal/cbioportal/archive/refs/tags/v5.4.10.tar.gz"
    extract_dir = os.path.join(tool_dir, "external_scripts")
    cbio_dir = os.path.join(extract_dir, "cbioportal-5.4.10")

    if not os.path.exists(cbio_dir):
        sys.stderr.write("Fetching and extracting cBioPortal validator scripts...\n")
        try:
            subprocess.run(f"mkdir -p {extract_dir} && curl -sL {repo_url} | tar xz -C {extract_dir}", shell=True, check=True)
            sys.stderr.write("cBioPortal repository successfully fetched and extracted.\n")
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"Error during fetching or extracting: {e}\n")
            raise
    
def setup_parser(tool_dir):
    parser = argparse.ArgumentParser(description="Run cBioPortal ETL pipeline")

    # General arguments
    parser.add_argument("--steps", nargs="+", type=str, help="Steps to execute (e.g., 1 2 3 or all)", choices=[str(i) for i in range(1, 6)] + ["all"])
    # Arguments for Step 1 - get_study_metadata.py
    parser.add_argument("-db", "--db-ini", required=True, help="Database config file")
    parser.add_argument("-p", "--profile", default="postgresql", help="Profile name (default: postgresql)")
    parser.add_argument("-mc", "--meta-config", required=False, help="Metadata configuration file. Default: value inputted for --study + '_case_meta_config.json'")
    parser.add_argument("-r", "--ref-dir", default=os.path.join(tool_dir, "REFS"), required=False, help="Reference directory. Defaults to tool's ref dir if not provided.")
    parser.add_argument("-a", "--all", required=False, action="store_true", help="Include all relevant files, not just status=active files, NOT RECOMMENDED")
    # Arguments for Step 2 - diff_studies.py
    parser.add_argument("-u", "--url", default="https://pedcbioportal.kidsfirstdrc.org/api/v2/api-docs", help="URL to search against")
    parser.add_argument("-s", "--study", required=True, help="Cancer study ID")
    parser.add_argument("-t", "--token", required=False, help="Token file obtained from Web API. Required if running Step 2")
    parser.add_argument("-ds", "--datasheet-sample", default=os.path.join(os.getcwd(), "datasheets/data_clinical_sample.txt"), help="File containing cBio-formatted sample metadata (default: datasheets/data_clinical_sample.txt from Step 1 output)")
    parser.add_argument("-dp", "--datasheet-patient", default=os.path.join(os.getcwd(), "datasheets/data_clinical_patient.txt"), help="File containing cBio-formatted patient metadata (default: datasheets/data_clinical_patient.txt from Step 1 output)")
    # Arguments for Step 3 - get_files_from_manifest.py
    parser.add_argument("-m", "--manifest", default="cbio_file_name_id.txt", help="Manifest file (default: cbio_file_name_id.txt from Step 1 output)")
    parser.add_argument("-f", "--file-types", required=False, help="Comma-separated file types to download")
    parser.add_argument("-at", "--aws-tbl", required=False, help="AWS table with bucket name and keys")
    parser.add_argument("-sp", "--sbg-profile", required=False, help="SBG profile name")
    parser.add_argument("-c", "--cbio", required=False, help="cBio manifest to limit downloads")
    parser.add_argument("-ao", "--active-only", required=False, action="store_true", help="Only include active files")
    parser.add_argument("-rm", "--rm-na", required=False, action="store_true", help="Remove entries where file_id and s3_path are NA")
    parser.add_argument("-d", "--debug", required=False, action="store_true", help="Enable debug mode")
    parser.add_argument("-o", "--overwrite", required=False, action="store_true", help="Overwrite files if they already exist")
    # Arguments for Step 4 - check_downloads.py
    parser.add_argument("-ms", "--manifest-subset", default="manifest_subset.tsv", required=False, help="Check that files were downloaded. Default: manifest_subset.tsv from Step 3")
    # Arguments for Step 5 - genomics_file_cbio_package_build.py
    parser.add_argument("-dc", "--data-config", required=False, help="Data processing configuration file. Default: value inputted for --study + '_data_processing_config.json'")
    parser.add_argument("-dgd", "--dgd-status", default="kf", required=False, choices=["both", "kf", "dgd"], help="Flag to determine load will have pbta/kf + dgd(both), kf/pbta only(kf), dgd-only(dgd). Default: kf")
    parser.add_argument("-l", "--legacy", required=False, action="store_true", help="Enable legacy mode")

    return parser


def main():
    tool_dir = os.path.dirname(os.path.abspath(__file__))
    configs_dir = os.path.join(tool_dir, "STUDY_CONFIGS")
    
    parser = setup_parser(tool_dir)
    args = parser.parse_args()
    args.meta_config = args.meta_config or os.path.join(configs_dir, f"{args.study}_case_meta_config.json")
    args.data_config = args.data_config or os.path.join(configs_dir, f"{args.study}_data_processing_config.json")

    steps = {
        "1": lambda: get_study_metadata(args),
        "2": lambda: diff_studies(args),
        "3": lambda: get_files_from_manifest(args),
        "4": lambda: check_downloads(args),
        "5": lambda: genomics_file_cbio_package_build(args),
    }

    selected_steps = args.steps
    if "all" in selected_steps: selected_steps = list(steps.keys())
    if "5" in selected_steps: fetch_validator_scripts(tool_dir)

    for step_number in selected_steps:
        if step_number in steps:
            print(f"\nRunning Step {step_number}...")
            steps[step_number]()
            print(f"Step {step_number} completed successfully.\n")
        else:
            print(f"Error in step {step_number}.")


if __name__ == "__main__":
    main()