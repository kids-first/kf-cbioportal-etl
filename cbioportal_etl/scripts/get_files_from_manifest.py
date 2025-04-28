#!/usr/bin/env python3
"""Download using a manifest files from AWS or SBG Platform."""

import argparse
import concurrent.futures
import os
import sys
from time import sleep
from typing import TYPE_CHECKING, Any

import boto3
import botocore
import pandas as pd
import sevenbridges as sbg
import urllib3
from sevenbridges.errors import SbgError
from sevenbridges.http.error_handlers import maintenance_sleeper, rate_limit_sleeper

if TYPE_CHECKING:
    from numpy import ndarray


def download_aws(
    file_type: str, key_dict: dict[str, dict[str, Any]], overwrite: bool, err_types: dict[str, int]
) -> None:
    """Download from AWS to file_type dir if AWS table was given.

    Args:
        file_type: String representation of genomic file type from ETL file
        key_dict: Dict with AWS key, profile, and s3 paths if applicable
        overwrite: Flag to overwrite existing files
        err_types: Dict to track any file download errors
    """
    for key in key_dict:
        current = key_dict[key]["manifest"]
        dl_client = key_dict[key]["dl_client"]
        sub_file_list = list(current.loc[current["file_type"] == file_type, "s3_path"])
        print(
            f"Grabbing {len(sub_file_list)} files using key {key}",
            file=sys.stderr,
        )
        for loc in sub_file_list:
            out: str = f"{file_type}/{loc.split('/')[-1]}"
            if overwrite or not os.path.isfile(out):
                parse_url = urllib3.util.parse_url(loc)
                try:
                    dl_client.download_file(
                        Bucket=parse_url.host,
                        Key=parse_url.path.lstrip("/"),
                        Filename=out,
                    )
                except Exception as e:
                    print(f"{e} could not download from {loc} using key {key}", file=sys.stderr)
                    sys.stderr.flush()
                    err_types["aws download"] += 1
            else:
                print(f"Skipping {out} it exists and overwrite not set", file=sys.stderr)
    print("Completed downloading files for " + file_type, file=sys.stderr)
    sys.stderr.flush()


def download_sbg(
    file_type: str,
    selected: pd.DataFrame,
    aws_tbl: str,
    api: sbg.Api,
    overwrite: bool,
    err_types: dict[str, int],
):
    """Download from SBG to file_type dir if SBG profile was given.

    Args:
        file_type: String representation of genomic file type from ETL file
        selected: Dataframe with filtered ETL entries
        api: SBG API object, if applicable
        overwrite: Flag to overwrite existing files
        err_types: Dict to track any file download errors
        aws_tbl: Table with AWS key info. Will trigger AWS downloads

    """
    sub_file_list: list = list(selected.loc[selected["file_type"] == file_type, "file_id"])
    # Sort of a "trust fall" that if aws bucket exists, skip SBG download
    if aws_tbl:
        sub_file_list = list(
            selected.loc[
                (selected["file_type"] == file_type) & (selected["s3_path"].isna()),
                "file_id",
            ]
        )
    for loc in sub_file_list:
        try:
            sbg_file: sbg.File = api.files.get(loc)
        except SbgError as e:
            print(
                f"Failed to get file with id {loc}. Will try once more in 3 seconds",
                file=sys.stderr,
            )
            print(e, file=sys.stderr)
            try:
                sleep(3)
                sbg_file = api.files.get(loc)
                print("Success on second try!", file=sys.stderr)
            except SbgError as e:
                print("Failed on second attempt", file=sys.stderr)
                err_types["sbg get"] += 1
        out = f"{file_type}/{sbg_file.name}"
        if overwrite or not os.path.isfile(out):
            try:
                sbg_file.download(out)
            except Exception as e:
                err_types["sbg download"] += 1
                print(f"Failed to download file with id {loc}", file=sys.stderr)
                print(e, file=sys.stderr)
        else:
            print(f"Skipping {out} it exists and overwrite not set", file=sys.stderr)
    return 0


def mt_type_download(
    file_type: str,
    key_dict: dict[str, dict[str, Any]],
    selected: pd.DataFrame,
    api: sbg.Api | None,
    overwrite: bool,
    err_types: dict[str, int],
    sbg_profile: str,
    aws_tbl: str,
) -> None:
    """Download files from each desired file type at the same time.

    Picks the download protocol based on provided args
    Args:
        file_type: String representation of genomic file type from ETL file
        key_dict: Dict with AWS key, profile, and s3 paths if applicable
        selected: Dataframe with filtered ETL entries
        api: SBG API object, if applicable
        overwrite: Flag to overwrite existing files
        err_types: Dict to track any file download errors
        sbg_profile: Name of SBG credentials to use from implicit credentials file. Will trigger SBG downloads
        aws_tbl: Table with AWS key info. Will trigger AWS downloads

    """
    if len(selected.loc[selected["file_type"] == file_type]):
        print(f"Downloading {file_type} files", file=sys.stderr)
        try:
            os.makedirs(file_type, exist_ok=True)
        except Exception as e:
            print(
                "{} error while making directory for {}".format(e, file_type),
                file=sys.stderr,
            )
        if aws_tbl:
            download_aws(file_type, key_dict, overwrite, err_types)
        if sbg_profile:
            download_sbg(file_type, selected, aws_tbl, api, overwrite, err_types)
    else:
        print(
            "WARNING: No files of type {file_type} in which file_id and s3_path is not NA. Skipping!",
            file=sys.stderr,
        )
    sys.stderr.flush()


def run_py(args: argparse.Namespace) -> int:
    # concat multiple possible manifests
    print("Concatenating manifests", file=sys.stderr)
    sys.stderr.flush()
    manifest_list: list[str] = args.manifest.split(",")
    manifest_df_list: list[pd.DataFrame] = []
    for manifest in manifest_list:
        print(f"Processing {manifest}", file=sys.stderr)
        manifest_df_list.append(pd.read_csv(manifest, sep=None, na_values=["NA"]))
    # In the event that s3_path is empty, replace with str to trigger later sbg download
    manifest_concat: pd.DataFrame = pd.concat(manifest_df_list, ignore_index=True)
    # if using a cbio name file as manifest, drop conflicting column
    if "etl_file_type" in manifest_concat.columns:
        del manifest_concat["etl_file_type"]
    # change col names to lower case for input compatibility
    manifest_concat.columns = manifest_concat.columns.str.lower()
    # file_types is actually a requirement, so grab from table if not provided
    if hasattr(args, "file_types") and args.file_types is not None:
        file_types_list = args.file_types.split(",")
    else:
        print("No file types provided, using table values", file=sys.stderr)
        file_types_list = manifest_concat["file_type"].unique().tolist()
    # subset concatenated manifests
    print("Subsetting concatenated manifest", file=sys.stderr)
    sys.stderr.flush()
    selected: pd.DataFrame = manifest_concat[manifest_concat["file_type"].isin(file_types_list)]
    # if active only flag set, further subset
    if args.active_only:
        print("active only flag given, will limit to that", file=sys.stderr)
        selected = selected[selected["status"] == "active"]
    if args.rm_na:
        print(
            "Drop na flag given. Will drop rows in which file_id and s3_path are NA",
            file=sys.stderr,
        )
        selected = selected[~(selected["file_id"].isna() & selected["s3_path"].isna())]
    # if cbio manifest given, limit to that
    if args.cbio:
        print(
            "cBio manifest provided, limiting downloads to matching IDs",
            file=sys.stderr,
        )
        cbio_data = pd.read_csv(args.cbio, sep=None, na_values=["NA"])
        specimen_list: ndarray = cbio_data.affected_bs_id.unique()
        selected = selected[selected["biospecimen_id"].isin(specimen_list)]
    # remove vcfs as we only want mafs
    pattern_del: str = ".vcf.gz"
    file_filter: pd.Series[bool] = selected["file_name"].str.contains(pattern_del)
    selected = selected[~file_filter]
    # special exception when file name is openpedcan, nothing to download
    opc_filter: pd.Series[bool] = selected["file_name"].str.contains("openpedcan")
    selected = selected[~opc_filter]
    out_file: str = "manifest_subset.tsv"
    selected.to_csv(out_file, sep="\t", mode="w", index=False)
    if args.debug:
        print(
            f"Debug flag given. No downloads actually happen, just a manifest subset to preview",
            file=sys.stderr,
        )
        sys.exit(1)
    err_types: dict[str, int] = {"aws download": 0, "sbg get": 0, "sbg download": 0}
    # download files by type
    check: int = 0
    key_dict = None
    api = None
    if args.aws_tbl is not None:
        check = 1
        # from https://stackoverflow.com/questions/66041582/connection-pool-is-full-warning-while-reading-s3-objects-via-multiple-threads
        client_config = botocore.config.Config(max_pool_connections=128)
        # setting up a key dict that, for each aws key, has an associated sesion and manifest to download with
        key_dict: dict[str, dict[str, Any]] = {}
        bucket_errs = 0
        with open(args.aws_tbl) as kl:
            for line in kl:
                (bucket, key) = line.rstrip("\n").split("\t")
                if key not in key_dict:
                    key_dict[key] = {}
                    key_dict[key]["manifest"] = selected[selected["s3_path"].str.contains(bucket)]
                    key_dict[key]["session"] = boto3.Session(profile_name=key)
                    key_dict[key]["dl_client"] = key_dict[key]["session"].client(
                        "s3", config=client_config
                    )
                else:
                    key_dict[key]["manifest"] = pd.concat(
                        [
                            key_dict[key]["manifest"],
                            selected[selected["s3_path"].str.startswith(bucket)],
                        ],
                        ignore_index=True,
                    )
                # Test bucket access with that key, if it fails, print error then kill to not waste time
                parse_url = urllib3.util.parse_url(bucket)
                try:
                    key_dict[key]["dl_client"].list_objects(Bucket=parse_url.host)
                except Exception as e:
                    bucket_errs = 1
                    print(e, file=sys.stderr)
                    print(f"Bucket access ERROR: {bucket}\t{key}", file=sys.stderr)
        if bucket_errs:
            sys.exit(1)
    if args.sbg_profile is not None:
        check = 1
        config: sbg.Config = sbg.Config(profile=args.sbg_profile)
        api = sbg.Api(config=config, error_handlers=[rate_limit_sleeper, maintenance_sleeper])
    if not check:
        print("Please provide at least one of aws_tbl or sbg_profile", file=sys.stderr)
        sys.exit(1)

    with concurrent.futures.ThreadPoolExecutor(16) as executor:
        {
            executor.submit(
                mt_type_download,
                ftype,
                key_dict,
                selected,
                api,
                args.overwrite,
                err_types,
                args.sbg_profile,
                args.aws_tbl,
            ): ftype
            for ftype in file_types_list
        }
    flag: int = 0
    for protocol, count in err_types.items():
        if count:
            flag += count
            print(f"ERROR: {protocol} {count} failure(s) occurred", file=sys.stderr)
    if flag:
        sys.exit(1)
    else:
        print("Downloads complete. No errors caught", file=sys.stderr)
    return 0


def main():
    parser = argparse.ArgumentParser(description="Get all files for a project.")
    parser.add_argument(
        "-m",
        "--manifest",
        action="store",
        dest="manifest",
        help="csv list of of genomic file location manifests",
    )
    parser.add_argument(
        "-f",
        "--file-types",
        action="store",
        dest="fts",
        help="csv list of workflow types to download",
    )
    parser.add_argument(
        "-at",
        "--aws-tbl",
        action="store",
        dest="aws_tbl",
        help="Table with bucket name and keys to subset on",
    )
    parser.add_argument(
        "-sp",
        "--sbg-profile",
        action="store",
        dest="sbg_profile",
        help="sbg profile name. Leave blank if using AWS instead",
    )
    parser.add_argument(
        "-c",
        "--cbio",
        action="store",
        dest="cbio",
        help="Add cbio manifest to limit downloads. Do NOT use if using cbio_file_name_id file as the manifest",
    )
    parser.add_argument(
        "-ao",
        "--active-only",
        default=False,
        action="store_true",
        dest="active_only",
        help="Set to grab only active files. Recommended.",
    )
    parser.add_argument(
        "-rm",
        "--rm-na",
        default=False,
        action="store_true",
        dest="rm_na",
        help="Remove entries where file_id and s3_path are NA.",
    )
    parser.add_argument(
        "-d",
        "--debug",
        default=False,
        action="store_true",
        dest="debug",
        help="Just output manifest subset to see what would be grabbed",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        default=False,
        action="store_true",
        dest="overwrite",
        help="If set, overwrite if file exists",
    )

    args = parser.parse_args()
    run_py(args)


if __name__ == "__main__":
    main()
