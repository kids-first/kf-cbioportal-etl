#!/usr/bin/env python
"""Read in GATK and ControlFreeC CNV data and convert to cBio format."""

import argparse
import csv
import json
import os
import sys
from typing import IO

from get_file_metadata_helper import get_file_metadata
from pybedtools import BedTool

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths


def get_gene_cnv_dict(cnv_bed_obj: BedTool, ref_bed: BedTool, ploidy: int = 2, high_gain: int = 4) -> tuple[dict[str, int], dict[str, int]]:
        # Create staging object that has cnv chrom, cnv start, cnv end, cn
        # ref chr, ref start, ref end, ref gene size overlap
        annotated: BedTool = cnv_bed_obj.intersect(ref_bed, b=True, wo=True)
        # iterate though annotated bed obj as list of lists to store in dict with uniq gene: CN values
        raw_cnv_dict = {}
        cn_idx: int = 3
        gene_idx: int = 7
        size_idx: int = 8
        # high gain is intended to be a relative threshold above ploidy,
        # in which when CN is ABOVE that value, it's considered an amplification (GISTIC 2)
        # for example if high_gain is 4, and ploidy is 2, high_gain for this sample becomes 6.
        # If ploidy is 3, high_gain becomes 7.
        high_gain += ploidy
        for entry in annotated:
            gene, cn, size = entry[gene_idx], int(entry[cn_idx]), int(entry[size_idx])
            # For repeat entries, keep largest size with greatest non-neutral CN
            if gene not in raw_cnv_dict or (size > raw_cnv_dict[gene][0] and cn != ploidy) or (size == raw_cnv_dict[gene][0] and abs(ploidy - cn) > abs(ploidy - raw_cnv_dict[gene][1])):
                raw_cnv_dict[gene] = [size, cn]
        # Drop sizes from raw and used to populate GISTIC-style dic
        gistic_cnv_dict: dict[str, int] = {}
        for gene in raw_cnv_dict:
            raw_cnv_dict[gene] = cn = raw_cnv_dict[gene][1]
            if cn > ploidy:
                g_cn = 1 if cn <= high_gain else 2
            elif cn < ploidy:
                g_cn = -1 if cn > 0 else -2
            else:
                g_cn = 0
            gistic_cnv_dict[gene] = g_cn
        return raw_cnv_dict, gistic_cnv_dict


def read_and_process_ctrlfreec_seg(ctrlfreec_seg_fname: str, out_seg_io: IO, sample_id: str) -> None:
    with open(ctrlfreec_seg_fname) as f:
        cnv_reader: csv.reader._reader = csv.reader(f, delimiter="\t")
        _header = next(cnv_reader)
        for entry in cnv_reader:
             entry[0] = sample_id
             entry[1] = entry[1][3:]  # remove leading chr
             print(*[entry], sep = "\t", file=out_seg_io)


def get_ctrlfreec_ploidy(cfree_info_fname: str) -> int:
    """Process ControFreeC Info file for [loidy information.

    Used to set a default value for CN if sample does not have a one for a specific gene
    Args:
        cfree_info_fname: Filename of ControFreeC info file
    Returns: str representation of output ploidy
    """
    info: IO = open(cfree_info_fname)
    for datum in info:
        datum: str = datum.rstrip("\n")
        try:
            if datum.startswith("Output_Ploidy"):
                (_key, value) = datum.split("\t")
                return int(value)
        except Exception as e:
            print(f"WARN: {e} entry could not be split", file=sys.stderr)
    return 2


def read_and_process_ctrlfreec_pval(pval_fname: str, ref_bed: BedTool, ploidy: int = 2) -> None:
    PVALUE: float = 0.05
    with open(pval_fname) as f:
        cnv_reader: csv.reader._reader = csv.reader(f, delimiter="\t")
        header = next(cnv_reader)
        wilcox_idx: int = header.index("WilcoxonRankSumTestPvalue")
        ks_idx: int = header.index("KolmogorovSmirnovPvalue")
        # create bed str representation of CNV data, filtering by pvalue
        cnv_filtered_as_bed: str = "\n".join(["\t".join(bed[0:4]) for bed in cnv_reader if float(bed[wilcox_idx]) <= PVALUE and float(bed[ks_idx]) <= PVALUE])
        cnv_bed_obj: BedTool = BedTool(cnv_filtered_as_bed, from_string=True)
        raw_cnv_dict, gistic_cnv_dict = get_gene_cnv_dict(cnv_bed_obj, ref_bed, ploidy)


def read_and_process_gatk_cnv(gatk_seg_fname: str, ref_bed: BedTool, sample_id: str, out_seg_io: IO):
    with open(gatk_seg_fname) as f:
        cnv_reader: csv.reader._reader = csv.reader(f, delimiter="\t")
        # Need to skip ridiculous interval list-style headers
        cnv_list = [entry for entry in cnv_reader if not entry[0].startswith("@")]
        # Drop chrM, convert MEAN_LOG2_COPY_RATIO to CN and rm leading chr
        # Also output to merged seg file while in the process
        cnv_as_bed = ""
        for entry in cnv_list[1:]:
            if entry[0] != "chrM":
                entry[0] = entry[0][3:]
                print(*[sample_id, *entry[0:4]], sep="\t", file=out_seg_io)
                # ratio is typically log2(cn/ploidy)
                entry[-2] = round(pow(2, float(entry[-2])) * 2)
                cnv_as_bed += f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[-2]}\n"
        cnv_bed_obj: BedTool = BedTool(cnv_as_bed, from_string=True)
        raw_cnv_dict, gistic_cnv_dict = get_gene_cnv_dict(cnv_bed_obj, ref_bed)

def main():
    parser = argparse.ArgumentParser(
        description="Merges files in gene <tab> entrez id <tab> copy number format into a genes-by-sample copy number table"
    )
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        dest="table",
        help="Table with cbio project, kf bs ids, cbio IDs, and file names",
    )

    parser.add_argument(
        "-j",
        "--config",
        action="store",
        dest="config_file",
        help="json config file with data types and data locations",
    )


    args = parser.parse_args()
    TOOL_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    with open(args.config_file) as f:
        config_data = json.load(f)
    config_data: dict = resolve_config_paths(config_data, TOOL_DIR)
    # for metadata dicts, cbio_project is primary key,  cbio_sample_name is secondary key.
    # Tertiary keys are other attributes with string as values.

    cnv_meta: dict[str, dict[str, dict[str, str]]] = get_file_metadata(args.table, "cnv")
    info_meta: dict[str, dict[str, dict[str, str]]] = get_file_metadata(args.table, "info")
    seg_meta: dict[str, dict[str, dict[str, str]]] = get_file_metadata(args.table, "seg")

    ref_bed = BedTool(config_data["bed_genes"])
    out_dir: str = "merged_cnvs"
    os.makedirs(out_dir, exist_ok=True)
    # cnv meta should have aLL IDs to be processed
    for project, samp_data in cnv_meta.items():
        # initialize out_seg_io here since it can be appended on-the-fly
        out_seg_io: IO = open(f"{out_dir}/{project}.merged_seg.txt", "w")
        print("ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean", file=out_seg_io)
        for cbio_sample, cbio_samp_attr in samp_data.items():
            if cbio_samp_attr["manifest_type"] == "ctrlfreec_pval":
                ploidy = 2
                if project in info_meta and cbio_sample in info_meta[project]:
                    ploidy = get_ctrlfreec_ploidy(info_meta[project][cbio_sample]["fname"])
                else:
                    print(f"WARNING: No info file for {cbio_sample} in {project}, using default ploidy of 2", file=sys.stderr)
                read_and_process_ctrlfreec_seg(ctrlfreec_seg_fname=seg_meta[project][cbio_sample]["fname"], out_seg_io=out_seg_io, sample_id=cbio_sample)
                read_and_process_ctrlfreec_pval(pval_fname=cbio_samp_attr["fname"], ref_bed=ref_bed, ploidy=ploidy)
            else:
                read_and_process_gatk_cnv(gatk_seg_fname=cbio_samp_attr["fname"], ref_bed=ref_bed, sample_id=cbio_sample, out_seg_io=out_seg_io)
        out_seg_io.close()

if __name__ == "__main__":
    main()
