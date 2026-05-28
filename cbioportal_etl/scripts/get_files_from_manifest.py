#!/usr/bin/env python3
"""Download using a manifest files from SBG Platform."""

import argparse
import concurrent.futures
import os
import sys
from math import ceil
from typing import TYPE_CHECKING, Any

import pandas as pd
import sevenbridges as sbg
from sevenbridges.errors import SbgError
from sevenbridges.http.error_handlers import maintenance_sleeper, rate_limit_sleeper

from cbioportal_etl.scripts.url_download_helper import parallel_download, small_download

if TYPE_CHECKING:
    from numpy import ndarray


def sbg_download_with_retry(file_obj: sbg.File, out: str, retries: int = 5, delay: int=3) -> None:
    """Download a file from SBG with retry logic.

    Will attempt to download a file up to `retries` times, with an exponential
    backoff delay between attempts.

    Args:
        file_obj: SBG File object to download
        out: Output file path
        retries: Number of retry attempts
        delay: Delay between retries in seconds

    """
    try:
        file_url = file_obj.download_info().url
        chunk_size: int = 32 * 1024 * 1024
        total_size: int = file_obj.size
        if total_size <= chunk_size:
            small_download(file_url, out, retries, delay)
        else:
            # set max worker size to ensure it doesn't exceed total file size for range requests
            max_workers: int = min(12, ceil(total_size / chunk_size))
            parallel_download(file_url, out, total_size, num_workers=max_workers,
                              chunk_size=chunk_size)

    except (Exception, SbgError) as e:
        print(f"Failed to download {out} after {retries} attempts from url: {file_url}",
              file=sys.stderr)
        dl_err_msg = f"Download failed for {file_obj.id} to location {out} with error {e}"
        raise Exception(dl_err_msg) from e


def download_sbg(
    file_type: str,
    selected: pd.DataFrame,
    api: sbg.Api,
    overwrite: bool,
    err_types: dict[str, int],
) -> dict[str, int]:
    """Download from SBG to file_type dir if SBG profile was given.

    Args:
        file_type: String representation of genomic file type from ETL file
        selected: Dataframe with filtered ETL entries
        api: SBG API object, if applicable
        overwrite: Flag to overwrite existing files
        err_types: Dict to track any file download errors

    """
    # get file id file name pairs from manifest
    sub_df = selected.loc[selected["file_type"] == file_type, ["file_id", "file_name", "s3_path"]]
    total_files = len(sub_df)

    batch_size = 100
    for batch_start_idx in range(0, len(sub_df), batch_size):
        print(f"Processed {file_type} files {batch_start_idx} out of {total_files}", file=sys.stderr)
        batch = sub_df.iloc[batch_start_idx : batch_start_idx + batch_size]
        batch_ids = batch["file_id"].tolist()
        batch_names = batch["file_name"].tolist()
        try:
            bulk_files = api.files.bulk_get(batch_ids)
            for j, file_obj in enumerate(bulk_files):
                e_type = "sbg get"
                if file_obj.valid:
                    out = f"{file_type}/{batch_names[j]}"
                    if not os.path.isfile(out) or overwrite:
                        e_type = "sbg download"
                        sbg_download_with_retry(file_obj.resource, out)
                    else:
                        print(f"Skipping {out} it exists and overwrite not set", file=sys.stderr)
                else:
                    print(f"File ID {batch_ids[j]} is not valid. Skipping download.",
                          file=sys.stderr)
                    err_types[e_type] += 1
        except sbg.SbgError as e:
            if e_type == "sbg download":
                print(f"Download failed for file {batch_names[j]}: {e}", file=sys.stderr)
            else:
                print(f"Bulk get failed for batch starting at index {batch_start_idx}: {e}",
                      file=sys.stderr)
            err_types[e_type] += 1
        except Exception as e:
            print(f"Unexpected error for batch starting at index {batch_start_idx}: {e}",
                  file=sys.stderr)
    print("Completed downloading files for " + file_type, file=sys.stderr)
    return err_types


def mt_type_download(
    file_type: str,
    selected: pd.DataFrame,
    api: sbg.Api | None,
    overwrite: bool,
    err_types: dict[str, int],
    sbg_profile: str,
) -> dict[str, int]:
    """Download files from each desired file type at the same time.

    Picks the download protocol based on provided args
    Args:
        file_type: String representation of genomic file type from ETL file
        selected: Dataframe with filtered ETL entries
        api: SBG API object, if applicable
        overwrite: Flag to overwrite existing files
        err_types: Dict to track any file download errors
        sbg_profile: Name of SBG credentials to use from implicit credentials file. Will trigger SBG downloads

    """
    if len(selected.loc[selected["file_type"] == file_type]):
        print(f"Downloading {file_type} files", file=sys.stderr)
        try:
            os.makedirs(file_type, exist_ok=True)
        except Exception as e:
            print(
                f"{e} error while making directory for {file_type}",
                file=sys.stderr,
            )
        if sbg_profile:
            err_types.update(download_sbg(file_type, selected, api, overwrite, err_types))
    else:
        print(
            f"WARNING: No files of type {file_type} in which file_id and s3_path is not NA. Skipping!",
            file=sys.stderr,
        )
    sys.stderr.flush()
    return err_types


def flag_adjustments(args: argparse.Namespace, selected: pd.DataFrame) -> pd.DataFrame:
    """Make adjustments to the manifest subset based on provided flags."""
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
        cbio_data = pd.read_csv(args.cbio, sep="\t", na_values=["NA"])
        specimen_list: ndarray = cbio_data.affected_bs_id.unique()
        selected = selected[selected["biospecimen_id"].isin(specimen_list)]
    return selected


def run_py(args: argparse.Namespace) -> int:
    """Run the main logic of the script."""
    # concat multiple possible manifests
    print("Concatenating manifests", file=sys.stderr)
    sys.stderr.flush()
    manifest_list: list[str] = args.manifest.split(",")
    manifest_df_list: list[pd.DataFrame] = []
    for manifest in manifest_list:
        print(f"Processing {manifest}", file=sys.stderr)
        manifest_df_list.append(pd.read_csv(manifest, sep="\t", na_values=["NA"]))
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
    selected = flag_adjustments(args, selected)
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
            "Debug flag given. No downloads actually happen, just a manifest subset to preview",
            file=sys.stderr,
        )
        sys.exit(1)
    err_types: dict[str, int] = {"sbg get": 0, "sbg download": 0}
    # download files by type
    if args.sbg_profile is not None:
        config: sbg.Config = sbg.Config(profile=args.sbg_profile)
        api = sbg.Api(config=config, error_handlers=[rate_limit_sleeper, maintenance_sleeper])
    else:
        print("Please provide sbg_profile", file=sys.stderr)
        sys.exit(1)

    with concurrent.futures.ThreadPoolExecutor(16) as executor:
        futures ={
            executor.submit(
                mt_type_download,
                ftype,
                selected,
                api,
                args.overwrite,
                err_types,
                args.sbg_profile,
            ): ftype
            for ftype in file_types_list
        }
        for fut in concurrent.futures.as_completed(futures):
            err_types.update(fut.result())
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


def main() -> None:
    """Parse args and run script."""
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
        help="csv list of workflow types to download",
    )
    parser.add_argument(
        "-sp",
        "--sbg-profile",
        action="store",
        dest="sbg_profile",
        help="sbg profile name",
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
