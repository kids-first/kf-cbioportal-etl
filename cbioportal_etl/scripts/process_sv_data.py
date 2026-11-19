#!/usr/bin/env python3
"""Converts RNA fusion, DNA annotSV and/or rMATS splice inputs into cBio SV files.

Any combination of the three sources can be run in one invocation; rows from all
selected sources are appended to the same per-project <project>.svs.txt.
"""

import argparse
import logging

from annotsv_to_cbio_sv import COLUMNS as DNA_COLUMNS
from annotsv_to_cbio_sv import iter_annotsv_chunks
from convert_fusion_as_sv import COLUMNS as FUSION_COLUMNS
from convert_fusion_as_sv import fusion_setup_outdir_metadata, iter_fusion_chunks
from rmats_splice_to_cbio_sv import COLUMNS as SPLICE_COLUMNS
from rmats_splice_to_cbio_sv import DEFAULT_COMMENT, iter_rmats_chunks
from sv_common import SVWriter, configure_logging, sv_setup_outdir_metadata


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser."""
    parser = argparse.ArgumentParser(
        description="Convert fusion, annotSV and/or rMATS splice results to cbio SV format."
    )
    parser.add_argument("-t", "--table", required=True,
                        help="Table with cbio project, kf bs ids, cbio IDs, and file names")
    parser.add_argument("-o", "--out-dir", default="merged_sv",
                        help="Result output dir. Default is merged_sv")
    parser.add_argument("-a", "--append", action=argparse.BooleanOptionalAction,
                        help="Append to existing project files, reusing their header")
    parser.add_argument("-s", "--skip-symlink", action=argparse.BooleanOptionalAction,
                        help="Skip sym links if repeating the process with the same input files")
    # RNA fusion
    parser.add_argument("-f", "--fusion-results",
                        help="annoFuse results dir OR openX merged fusion file")
    parser.add_argument("-m", "--mode", choices=["openX", "kfprod", "dgd"],
                        help="Fusion source (required with --fusion-results)")
    # DNA SV
    parser.add_argument("--process-sv-results", action=argparse.BooleanOptionalAction,
                        help="Process DNA SV results from annotSV listed in the table")
    parser.add_argument("-l", "--link-input-dir", default="annotSV_results/",
                        help="DNA SV directory name to store annotSV symlinks")
    # rMATS splice
    parser.add_argument("--process-splice-results", action=argparse.BooleanOptionalAction,
                        help="Process rMATS splice results listed in the table")
    parser.add_argument("--splice-link-input-dir", default="rmats_results/",
                        help="rMATS splice directory name to store symlinks")
    parser.add_argument("--comment", default=DEFAULT_COMMENT,
                        help="rMATS and GENCODE version information for the splice Comments column")
    return parser


def main() -> None:
    """Convert any combination of fusion, annotSV and rMATS results to cBio SV format."""
    parser = build_parser()
    args = parser.parse_args()
    if not (args.fusion_results or args.process_sv_results or args.process_splice_results):
        parser.error("Specify at least one of --fusion-results, --process-sv-results, --process-splice-results.")
    if args.fusion_results and not args.mode:
        parser.error("Need to set -m mode when processing fusion data. Opts are kfprod, openX, dgd.")
    configure_logging()
    logger = logging.getLogger(__name__)

    header: list[str] = []
    if args.process_sv_results:
        header += DNA_COLUMNS
    if args.fusion_results:
        header += FUSION_COLUMNS
    if args.process_splice_results:
        header += SPLICE_COLUMNS
    writer = SVWriter(args.out_dir, list(dict.fromkeys(header)), append=bool(args.append))

    projects: set[str] = set()
    if args.process_sv_results:
        logger.info("DNA SV flag given, processing")
        meta = sv_setup_outdir_metadata(args.link_input_dir, args.out_dir, args.table, "sv", args.skip_symlink)
        projects.update(meta.cbio_project.unique())
        for project, df in iter_annotsv_chunks(args.link_input_dir, meta):
            writer.write(project, df)
    if args.fusion_results:
        logger.info("RNA fusion data given, processing")
        meta = fusion_setup_outdir_metadata(args.mode, args.out_dir, args.table)
        projects.update(meta.cbio_project.unique())
        for project, df in iter_fusion_chunks(args.fusion_results, args.mode, meta):
            writer.write(project, df)
    if args.process_splice_results:
        logger.info("rMATS splice data given, processing")
        meta = sv_setup_outdir_metadata(
            args.splice_link_input_dir, args.out_dir, args.table, "splice", args.skip_symlink
        )
        projects.update(meta.cbio_project.unique())
        for project, df in iter_rmats_chunks(args.splice_link_input_dir, meta, args.comment):
            writer.write(project, df)
    writer.finish(sorted(projects))
    logger.info("Done writing SV files for %d project(s) to %s", len(projects), args.out_dir)


if __name__ == "__main__":
    main()
