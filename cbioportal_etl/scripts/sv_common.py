#!/usr/bin/env python3
"""Shared helpers for the cBio structural variant (SV) converters.

Each converter exposes a generator yielding ``(cbio_project, formatted_df)`` chunks
(one input file, or one batch of samples, at a time). ``SVWriter`` appends those
chunks to ``<project>.svs.txt`` so memory use stays at roughly one chunk.
"""

import csv
import logging
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

NA_REP = ""  # how missing values (NaN or <NA>) are written


def configure_logging() -> None:
    """Configure root logging. Call from main() only, never at import time."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def sv_setup_outdir_metadata(
    link_input_dir: str, out_dir: str, table: str, ftype: str, skip: bool
) -> pd.DataFrame:
    """Create output dir and symlink dir, symlink inputs, and subset metadata by file type.

    Args:
        link_input_dir: Directory to store symlinks to the input results
        out_dir: Directory to store merged SV results
        table: Table with cbio project, kf bs ids, cbio IDs, and file names
        ftype: etl_file_type value to filter for (e.g. "sv", "splice")
        skip: Skip symlinking if repeating the process with the same input files
    Returns:
        Pandas dataframe with the metadata rows for ftype

    """
    logger = logging.getLogger(__name__)
    Path(link_input_dir).mkdir(exist_ok=True, parents=True)
    Path(out_dir).mkdir(exist_ok=True, parents=True)

    # ensure sample name is imported as str
    all_file_meta = pd.read_csv(table, sep="\t", dtype={"cbio_sample_name": str})
    subset = all_file_meta.loc[all_file_meta["etl_file_type"] == ftype]
    in_dirs: list[str] = [] if skip else subset.file_type.unique().tolist()
    logger.info(
        "%sSymlinking %s files from %s to %s",
        "Skipping " if skip else "", ftype, ",".join(subset.file_type.unique()), link_input_dir,
    )
    sym_errs = 0
    for dirname in in_dirs:
        if not Path(dirname).exists():
            err = f"Input dir {dirname} does not exist. Check if files were downloaded to the correct location"
            raise ValueError(err)
        for src in Path(dirname).resolve().iterdir():
            try:
                (Path(link_input_dir) / src.name).symlink_to(src)
            except OSError:
                logger.exception("Could not sym link %s in %s", src, dirname)
                sym_errs += 1
    # If symlink errors, stop here as data will be incomplete
    if sym_errs:
        err = f"Could not sym link {sym_errs} files, exiting!"
        raise FileNotFoundError(err)
    return subset


class SVWriter:
    """Append per-project chunks to ``<out_dir>/<project>.svs.txt``.

    The header is fixed up front, so chunks from different sources can be mixed.
    Chunks are reindexed to the header; columns a chunk lacks are written blank.
    In append mode the header is taken from the existing file instead.
    """

    def __init__(self, out_dir: str, header: list[str], *, append: bool = False) -> None:
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(exist_ok=True, parents=True)
        self.header = header
        self.append = append
        self._cols: dict[str, list[str]] = {}

    def _path(self, project: str) -> Path:
        return self.out_dir / f"{project}.svs.txt"

    def write(self, project: str, df: pd.DataFrame) -> None:
        """Write one chunk for a project, creating the file with a header on first use."""
        path = self._path(project)
        first = project not in self._cols
        if first:
            self._cols[project] = (
                list(pd.read_csv(path, sep="\t", nrows=0).columns) if self.append else self.header
            )
        write_header = first and not self.append
        df.reindex(columns=self._cols[project]).to_csv(
            path, sep="\t", index=False, quoting=csv.QUOTE_NONE, na_rep=NA_REP,
            mode="w" if write_header else "a", header=write_header,
        )

    def finish(self, projects: Iterable[str]) -> None:
        """Write header-only files for projects that never received any rows."""
        if self.append:
            return
        for project in projects:
            if project not in self._cols:
                self.write(project, pd.DataFrame(columns=self.header))