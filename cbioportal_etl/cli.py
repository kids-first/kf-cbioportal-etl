"""Collate args for cBio ETL."""

import argparse
import os
import sys

from cbioportal_etl.steps import run_etl


def main():
    parser = argparse.ArgumentParser(description="CBio ETL Command Line Tool")
    subparsers = parser.add_subparsers(dest="command", required=True)

    common_args = argparse.ArgumentParser(add_help=False)
    # Step 1 - generate_config.py
    common_args.add_argument(
        "-ctsv",
        "--config-tsv",
        required=False,
        help="Path to the input TSV file containing static values (preset)",
    )
    common_args.add_argument(
        "-stsv",
        "--study-tsv",
        required=False,
        help="Path to the input TSV file containing completed values for json generation",
    )
    common_args.add_argument("-s", "--study", required=True, help="Cancer study ID")
    # Step 2 - get_study_metadata.py
    common_args.add_argument("-db", "--db-ini", required=True, help="Database config file")
    common_args.add_argument(
        "-p", "--profile", default="postgresql", help="Profile name (default: postgresql)"
    )
    common_args.add_argument(
        "-sc",
        "--study-config",
        required=False,
        help="JSON file generated from Step 1 or path to custom JSON (automatically set when applicable)",
    )
    common_args.add_argument("-r", "--ref-dir", required=False, help="Reference directory")
    common_args.add_argument(
        "-a",
        "--all",
        required=False,
        action="store_true",
        help="Include all relevant files, not just status=active files",
    )
    # Step 4 - get_files_from_manifest.py
    common_args.add_argument(
        "-m",
        "--manifest",
        default="cbio_file_name_id.txt",
        help="Manifest file (automatically set when applicable)",
    )
    common_args.add_argument(
        "-f", "--file-types", required=False, help="Comma-separated file types to download"
    )
    common_args.add_argument(
        "-at", "--aws-tbl", required=False, help="AWS table with bucket name and keys"
    )
    common_args.add_argument("-sp", "--sbg-profile", required=True, help="SBG profile name")
    common_args.add_argument(
        "-c", "--cbio", required=False, help="cBio manifest to limit downloads"
    )
    common_args.add_argument(
        "-ao",
        "--active-only",
        required=False,
        action="store_true",
        help="Only include active files",
    )
    common_args.add_argument(
        "-rm",
        "--rm-na",
        required=False,
        action="store_true",
        help="Remove entries where file_id and s3_path are NA",
    )
    common_args.add_argument(
        "-d", "--debug", required=False, action="store_true", help="Enable debug mode"
    )
    common_args.add_argument(
        "-o",
        "--overwrite",
        required=False,
        action="store_true",
        help="Overwrite files if they already exist",
    )
    # Step 5 - check_downloads.py
    common_args.add_argument(
        "-ms",
        "--manifest-subset",
        default="manifest_subset.tsv",
        required=False,
        help="Check that files were downloaded (automatically set when applicable)",
    )
    # Step 6 - genomics_file_cbio_package_build.py
    common_args.add_argument(
        "-dgd", "--dgd-status", default="kf", choices=["both", "kf", "dgd"], help="Load options"
    )
    common_args.add_argument(
        "-ad",
        "--add-data",
        action="store_true",
        dest="add_data",
        help="Flag to indicate add_data mode (automatically set when applicable)",
    )
    common_args.add_argument(
        "-et", 
        "--expression-type", 
        action="store", 
        dest="expression_type", 
        choices=["TPM", "FPKM"], 
        default="TPM", 
        help="Which expression value to use: TPM or FPKM. Default is TPM."
    )
    common_args.add_argument(
        "-dmt", 
        "--default-match-type", 
        action="store", 
        dest="default_match_type", 
        choices=["polyA", "totalRNA", "none"], 
        default="none", 
        help="Default match type for samples with unknown RNA library type for z-score calculations. Use 'polyA' or 'totalRNA' to override fallback to intra-cohort z-score."
    ) 

    # Arguments exclusive to update (step 3 - diff_studies.py)
    update_args = argparse.ArgumentParser(add_help=False)
    update_args.add_argument(
        "-u",
        "--url",
        default="https://pedcbioportal.kidsfirstdrc.org/api/v2/api-docs",
        help="URL to search against",
    )
    update_args.add_argument(
        "-ds",
        "--datasheets",
        default="datasheets",
        required=False,
        help="Directory containing cBio-formatted metadata (automatically set when applicable)",
    )
    update_args.add_argument(
        "-t", "--token", required=True, help="Token file obtained from Web API"
    )

    # Import command (full study import)
    subparsers.add_parser(
        "import", parents=[common_args], help="Run import workflow (Steps 1, 2, 4, 5, 6)"
    )

    # Update command (incremental updates)
    subparsers.add_parser(
        "update",
        parents=[common_args, update_args],
        help="Run update workflow (Steps 1, 2, 3, 4, 5, 6)",
    )

    args: argparse.Namespace = parser.parse_args()

    if args.command == "import":
        run_etl(args, steps=["1", "2", "4", "5", "6"])
    elif args.command == "update":
        run_etl(args, steps=["1", "2", "3"])
        # If there are new patients to be added, run the rest of the ETL
        study_add_data_dir: str = f"{args.study}_add_data"
        cbio_manifest_file: str = f"{study_add_data_dir}/cbio_file_name_id.txt"
        if os.path.exists(study_add_data_dir) and os.path.isfile(cbio_manifest_file):
            print(f"\nDetected {study_add_data_dir} directory. Switching to it.\n")
            for path_arg in ["study_config", "file_types", "aws_tbl", "cbio"]:
                arg_val = getattr(args, path_arg, None)
                if arg_val and not os.path.isabs(arg_val):
                    setattr(args, path_arg, os.path.abspath(arg_val))

            os.chdir(study_add_data_dir)
            args.manifest = "cbio_file_name_id.txt"
            args.datasheets = "datasheets"
            args.add_data = True
            run_etl(args, steps=["4", "5", "6"])
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
