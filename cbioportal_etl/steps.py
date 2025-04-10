"""cBio ETL wrapper."""

import os
import subprocess
import sys

from cbioportal_etl.scripts.check_downloads import run_py as check_downloads
from cbioportal_etl.scripts.diff_studies import run_py as diff_studies
from cbioportal_etl.scripts.generate_config import run_py as generate_config
from cbioportal_etl.scripts.genomics_file_cbio_package_build import (
    run_py as genomics_file_cbio_package_build,
)
from cbioportal_etl.scripts.get_files_from_manifest import run_py as get_files_from_manifest
from cbioportal_etl.scripts.get_study_metadata import run_py as get_study_metadata


def fetch_validator_scripts(tool_dir: str) -> None:
    """Get validator scripts from MSKCC cbioportal repo."""
    repo_url: str = "https://github.com/cBioPortal/cbioportal/archive/refs/tags/v5.4.10.tar.gz"
    extract_dir: str = os.path.join(tool_dir, "external_scripts")
    cbio_dir: str = os.path.join(extract_dir, "cbioportal-5.4.10")

    if not os.path.exists(cbio_dir):
        print("Fetching and extracting cBioPortal validator scripts...", file=sys.stderr)
        try:
            subprocess.run(
                f"mkdir -p {extract_dir} && curl -sL {repo_url} | tar xz -C {extract_dir}",
                shell=True,
                check=True,
            )
            print("cBioPortal repository successfully fetched and extracted.", file=sys.stderr)
        except subprocess.CalledProcessError as e:
            print(f"Error during fetching or extracting: {e}", file=sys.stderr)
            raise


def fetch_dgd_fusion():
    """
    Download fusion-dgd.tsv.gz if it doesn't already exist in cwd
    """
    filename = "fusion-dgd.tsv.gz"
    url = "https://d3b-openaccess-us-east-1-prd-pbta.s3.us-east-1.amazonaws.com/open-targets/v15/fusion-dgd.tsv.gz"
    if not os.path.exists(filename):
        sys.stderr.write(f"{filename} not found. Downloading from AWS...\n")
        try:
            subprocess.run(
                f"curl -sSL -o {filename} {url}",
                shell=True,
                check=True
            )
            sys.stderr.write(f"Downloaded {filename} successfully.\n")
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"Failed to download {filename}: {e}\n")
            raise


def run_etl(args, steps):
    tool_dir = os.path.dirname(os.path.abspath(__file__))

    args.study_config = args.study_config or f"{args.study}_config.json"
    args.config_tsv = args.config_tsv or os.path.join(
        tool_dir, "STUDY_CONFIGS/all_studies_config_values.tsv"
    )
    args.ref_dir = args.ref_dir or os.path.join(tool_dir, "REFS")
    if args.dgd_status != "kf":
        fetch_dgd_fusion()

    steps_map = {
        "1": lambda: generate_config(args),
        "2": lambda: get_study_metadata(args),
        "3": lambda: diff_studies(args),
        "4": lambda: get_files_from_manifest(args),
        "5": lambda: check_downloads(args),
        "6": lambda: genomics_file_cbio_package_build(args),
    }

    if "6" in steps:
        fetch_validator_scripts(tool_dir)

    for step in steps:
        if step in steps_map:
            print(f"\nRunning Step {step}...")
            try:
                steps_map[step]()
                print(f"Step {step} completed successfully.\n")
            except Exception as e:
                print(f"Error in Step {step}: {e}", file=sys.stderr)
                sys.exit(1)
        else:
            print(f"Error: Invalid step {step}.")
            sys.exit(1)
