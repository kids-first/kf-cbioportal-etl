#!/usr/bin/env python3
"""Merge MAF Files.

Merge and filter MAF files to create an output with a singular header and
variants not normally suitable for cBio removed

"""

import argparse
import logging
import sys
from operator import itemgetter
from pathlib import Path
from typing import IO

import pandas as pd

logger = logging.getLogger(__name__)

MISSING = -1  # index of the "" sentinel appended to every row


def _column_index(header_line: str) -> dict[str, int]:
    """Map column name -> position for one MAF header line."""
    return {name: i for i, name in enumerate(header_line.rstrip("\n").split("\t"))}


def _make_output_picker(col: dict[str, int], print_header: list[str]):
    """Build a callable that pulls the output fields from a split row, in order.

    Columns absent from this file resolve to MISSING, which points at the
    trailing "" appended to each row.
    """
    indices = [col.get(name, MISSING) for name in print_header]
    if len(indices) == 1:  # itemgetter returns a bare item (not a tuple) for one index
        only = indices[0]
        return lambda row: (row[only],)
    return itemgetter(*indices)


def process_maf(
        maf_fn: Path,
        new_maf: IO,
        maf_exc: dict[str, int],
        tum_id: str,
        norm_id: str,
        print_header: list[str]) -> None:
    """Process a single MAF file, filtering out unwanted entries and writing to new MAF."""
    with maf_fn.open(buffering=1 << 20) as cur_maf:  # ~1MB reads
        next(cur_maf)  # first line (version/comment), skipped as in the original

        # --- one-time setup: all the readable structure lives here ---
        col = _column_index(next(cur_maf))
        tumor_idx = col["Tumor_Sample_Barcode"]
        norm_idx = col["Matched_Norm_Sample_Barcode"]
        class_idx = col["Variant_Classification"]
        gene_idx = col["Hugo_Symbol"]
        pick_output = _make_output_picker(col, print_header)

        # Filtering only needs the leading fields, so split just that far.
        n_leading_fields = max(class_idx, gene_idx) + 1
        write = new_maf.write

        # --- hot loop: kept flat on purpose (no per-row function calls) ---
        for line in cur_maf:
            head = line.split("\t", n_leading_fields)  # partial split: cheap reject
            variant_class = head[class_idx]
            is_excluded = variant_class in maf_exc and not (
                variant_class == "5'Flank" and head[gene_idx] == "TERT"
            )
            if is_excluded:
                continue

            row = line.rstrip("\n").split("\t")  # full split only for kept rows
            row[tumor_idx] = tum_id
            row[norm_idx] = norm_id
            row.append("")  # sentinel for MISSING columns
            write("\t".join(pick_output(row)) + "\n")


def process_tbl(study: str, file_meta_dict: list[list], print_head: str) -> None:
    """Process MAF file be project (cbio_dx).

    Unlikely that more than one project will be processed, pass each variant entry for filtering
    Args:
    cbio_dx: cBio project name
    file_meta_dict: Dict that has been subset by file type from ETL file
    print_head: Output header line for file output an dto guide output content
    """
    try:
        x = 0
        # project/disease name should be name of directory hosting datasheet
        logger.info("Processing %s project", study)
        with (out_dir / f"{study}.maf").open("w") as new_maf:
            new_maf.write(print_head)
            for fname, cbio_tum_id, cbio_norm_id in file_meta_dict:
                logger.info("Found relevant maf to process for %s %s %s", cbio_tum_id, cbio_norm_id, fname)
                process_maf(maf_dir / fname, new_maf, maf_exc, cbio_tum_id, cbio_norm_id, print_header)
                x += 1
            logger.info("Completed processing %d entries in %s", x, study)
    except Exception:
        logger.exception("Failed while processing %s", study)
        sys.exit(1)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    parser = argparse.ArgumentParser(
        description="Merge and filter MAFs using ETL table withe file locations.",
    )
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        dest="table",
        help="Table with cbio project, kf bs ids, cbio IDs, and file names",
    )
    parser.add_argument(
        "-i",
        "--header",
        action="store",
        dest="header",
        help="File with maf header only",
    )
    parser.add_argument(
        "-s",
        "--skip-symlink",
        action="store_true",
        help="Flag to skip sym links if repeating the process with the same input files",
    )


    args = parser.parse_args()
    # Create symlinks to mafs in one place for ease of processing
    maf_dir: Path = Path("MAFS")
    all_file_meta = pd.read_csv(args.table, sep="\t", dtype={"cbio_sample_name": str})
    maf_subset = all_file_meta[all_file_meta["etl_file_type"] == "maf"].copy()
    study = maf_subset["cbio_project"].iloc[0]
    maf_list = maf_subset[["file_name", "cbio_sample_name", "cbio_matched_normal_name"]].fillna("").astype(str).drop_duplicates().values.tolist()
    if args.skip_symlink:
        logger.info("Skipping symlinking of MAFs to %s", maf_dir)
    else:
        maf_dirs_in: list[str] = maf_subset.file_type.unique().tolist()
        logger.info("Symlinking maf files from %s to %s", ",".join(maf_dirs_in), maf_dir)
        maf_dir.mkdir(exist_ok=True)
        sym_errs: int = 0
        for dirname in maf_dirs_in:
            if not Path(dirname).exists():
                logger.error(
                    "Input MAF dir %s does not exist. Check if files were downloaded to the correct location",
                    dirname,
                )
                sys.exit(2)
            for src in Path(dirname).resolve().iterdir():
                try:
                    (maf_dir / src.name).symlink_to(src)
                except OSError:
                    logger.exception("Could not sym link %s in %s", src.name, dirname)
                    sym_errs += 1
        # If symlink errors, stop here as data will be incomplete
        if sym_errs:
            logger.error("Could not sym link %d files, exiting!", sym_errs)
            sys.exit(2)

    with Path(args.header).open() as head_fh:
        print_head: str = next(head_fh)
        cur_head: str = next(head_fh)
        print_header: list[str] = cur_head.rstrip("\n").split("\t")
        eid_idx: int = print_header.index("Entrez_Gene_Id")
        print_header.pop(eid_idx)

    print_head += "\t".join(print_header) + "\n"
    maf_exc: dict[str, int] = {
        "Silent": 0,
        "Intron": 0,
        "IGR": 0,
        "3'UTR": 0,
        "5'UTR": 0,
        "3'Flank": 0,
        "5'Flank": 0,
        "RNA": 0,
    }
    out_dir: Path = Path("merged_mafs")
    out_dir.mkdir(exist_ok=True, parents=True)
    # iterating through projects that are the first key in dict
    process_tbl(study, maf_list, print_head)

    logger.info("Done processing MAFs for %s", study)
