#!/usr/bin/env python3
"""Download using a manifest files from SBG Platform."""

import argparse
import concurrent.futures
import logging
import os
import sys
from typing import TYPE_CHECKING

import pandas as pd
import sevenbridges as sbg
from sevenbridges.errors import SbgError
from sevenbridges.http.error_handlers import maintenance_sleeper, rate_limit_sleeper

from cbioportal_etl.scripts.url_download_helper import parallel_download, small_download

if TYPE_CHECKING:
    from numpy import ndarray
from threading import Lock


class DownloadTracker:
    """Class to track download status of files across threads."""

    def __init__(self) -> None:
        """Initialize the tracker with empty lists and a lock."""
        self.failed = []
        self.invalid = []
        self.success = []
        self.lock = Lock()

    def add_failed(self, file_id: str, path: str, error: Exception) -> None:
        """Add a failed download entry to the tracker."""
        with self.lock:
            self.failed.append((file_id, path, str(error)))

    def add_invalid(self, file_id: str) -> None:
        """Add an invalid file ID to the tracker."""
        with self.lock:
            self.invalid.append(file_id)

    def add_success(self, file_id: str, path: str) -> None:
        """Add a successful download entry to the tracker."""
        with self.lock:
            self.success.append((file_id, path))

def sbg_download_with_retry(
        file_obj: sbg.File,
        out: str,
        tracker: DownloadTracker,
        retries: int = 5,
        delay: int=3,
        ) -> None:
    """Download a file from SBG with retry logic.

    Will attempt to download a file up to `retries` times, with an exponential
    backoff delay between attempts.

    Args:
        file_obj: SBG File object to download
        out: Output file path
        tracker: DownloadTracker object to track download status
        retries: Number of retry attempts
        delay: Delay between retries in seconds

    """
    logger = logging.getLogger(__name__)
    file_id = file_obj.id
    try:
        file_url = file_obj.download_info().url
        chunk_size: int = 32 * 1024 * 1024
        total_size: int = file_obj.size
        if total_size <= chunk_size:
            small_download(file_url, out, retries, delay)
        else:
            parallel_download(file_url, out, total_size, num_workers=12,
                              chunk_size=chunk_size)
        tracker.add_success(file_id, out)
    except (Exception, SbgError) as e:
        logger.exception("Failed to download %s after %d attempts from url: %s", out, retries, file_url)
        tracker.add_failed(file_id, out, e)
        dl_err_msg = f"Download failed for {file_obj.id} to location {out} with error {e}"
        raise Exception(dl_err_msg) from e


def download_sbg(
    file_type: str,
    selected: pd.DataFrame,
    api: sbg.Api,
    overwrite: bool,
    tracker: DownloadTracker,
) -> None:
    """Download from SBG to file_type dir if SBG profile was given.

    Args:
        file_type: String representation of genomic file type from ETL file
        selected: Dataframe with filtered ETL entries
        api: SBG API object, if applicable
        overwrite: Flag to overwrite existing files
        tracker: DownloadTracker object to track download status across threads

    """
    logger = logging.getLogger(__name__)
    # get file id file name pairs from manifest
    sub_df = selected.loc[selected["file_type"] == file_type, ["file_id", "file_name", "s3_path"]]
    total_files = len(sub_df)

    batch_size = 100
    for batch_start_idx in range(0, len(sub_df), batch_size):
        logger.info("Processed %s files %d out of %d", file_type, batch_start_idx, total_files)
        batch = sub_df.iloc[batch_start_idx : batch_start_idx + batch_size]
        batch_ids = batch["file_id"].tolist()
        batch_names = batch["file_name"].tolist()
        try:
            bulk_files = api.files.bulk_get(batch_ids)
            for j, file_obj in enumerate(bulk_files):
                if file_obj.valid:
                    out = f"{file_type}/{batch_names[j]}"
                    if not os.path.isfile(out) or overwrite:
                        sbg_download_with_retry(file_obj.resource, out, tracker)

                    else:
                        logger.info("Skipping %s it exists and overwrite not set", out)
                else:
                    logger.warning("File ID %s is not valid. Skipping download.", batch_ids[j])
                    tracker.add_invalid(batch_ids[j])
        except Exception as e:
            logger.exception("Unexpected error for batch starting at index %d: %s", batch_start_idx, e)
    logger.info("Completed downloading files for %s", file_type)


def mt_type_download(
    file_type: str,
    selected: pd.DataFrame,
    api: sbg.Api,
    overwrite: bool,
    tracker: DownloadTracker,
) -> None:
    """Download files from each desired file type at the same time.

    Picks the download protocol based on provided args
    Args:
        file_type: String representation of genomic file type from ETL file
        selected: Dataframe with filtered ETL entries
        api: SBG API object, if applicable
        overwrite: Flag to overwrite existing files
        tracker: DownloadTracker object to track download status across threads

    """
    logger = logging.getLogger(__name__)
    if len(selected.loc[selected["file_type"] == file_type]):
        logger.info("Downloading %s files", file_type)
        try:
            os.makedirs(file_type, exist_ok=True)
            download_sbg(file_type, selected, api, overwrite, tracker)
        except Exception as e:
            logger.exception("error while making directory for %s", file_type)
    else:
        logger.warning(
            "No files of type %s in which file_id and s3_path is not NA. Skipping!", file_type
        )
    sys.stderr.flush()


def flag_adjustments(args: argparse.Namespace, selected: pd.DataFrame) -> pd.DataFrame:
    """Make adjustments to the manifest subset based on provided flags."""
    # if active only flag set, further subset
    logger = logging.getLogger(__name__)
    if args.active_only:
        logger.info("active only flag given, will limit to that")
        selected = selected[selected["status"] == "active"]
    if args.rm_na:
        logger.info("Drop na flag given. Will drop rows in which file_id and s3_path are NA")
        selected = selected[~(selected["file_id"].isna() & selected["s3_path"].isna())]
    # if cbio manifest given, limit to that
    if args.cbio:
        logger.info("cBio manifest provided, limiting downloads to matching IDs")
        cbio_data = pd.read_csv(args.cbio, sep="\t", na_values=["NA"])
        specimen_list: ndarray = cbio_data.affected_bs_id.unique()
        selected = selected[selected["biospecimen_id"].isin(specimen_list)]
    return selected


def run_py(args: argparse.Namespace) -> int:
    """Run the main logic of the script."""
    # concat multiple possible manifests
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        )
    logger = logging.getLogger(__name__)
    logger.info("Concatenating manifests")
    manifest_list: list[str] = args.manifest.split(",")
    manifest_df_list: list[pd.DataFrame] = []
    for manifest in manifest_list:
        logger.info("Processing %s", manifest)
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
        logger.info("No file types provided, using table values")
        file_types_list = manifest_concat["file_type"].unique().tolist()
    # subset concatenated manifests
    logger.info("Subsetting concatenated manifest")
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
        logger.info("Debug flag given. No downloads actually happen, just a manifest subset to preview")
        sys.exit(0)
    # download files by type
    if args.sbg_profile is None:
        logger.error("Please provide sbg_profile")
        sys.exit(1)
    config: sbg.Config = sbg.Config(profile=args.sbg_profile)
    api = sbg.Api(config=config, error_handlers=[rate_limit_sleeper, maintenance_sleeper])

    tracker = DownloadTracker()
    with concurrent.futures.ThreadPoolExecutor(16) as executor:
        futures ={
            executor.submit(
                mt_type_download,
                ftype,
                selected,
                api,
                args.overwrite,
                tracker,
            ): ftype
            for ftype in file_types_list
        }

        for future in concurrent.futures.as_completed(futures):
            future.result()
    logger.info("======== DOWNLOAD SUMMARY ========")
    logger.info("Successful downloads: %d", len(tracker.success))
    logger.info("Failed downloads: %d", len(tracker.failed))
    logger.info("Invalid file IDs: %d", len(tracker.invalid))

    if tracker.failed:
        logger.info("---- FAILED DOWNLOADS ----")
        for file_id, path, err in tracker.failed:
            logger.info("%s -> %s | %s", file_id, path, err)

    if tracker.invalid:
        logger.info("---- INVALID FILE IDs ----")
        for file_id in tracker.invalid:
            logger.info("%s", file_id)

    if not tracker.failed and not tracker.invalid:
        logger.info("All files downloaded successfully with valid IDs!")
        return 0
    logger.error("Some files failed to download or had invalid IDs. Please review the logs for details.")
    return 1


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
