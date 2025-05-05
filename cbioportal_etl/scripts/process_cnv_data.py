#!/usr/bin/env python
"""Read in GATK and ControlFreeC CNV data and convert to cBio format."""

import argparse
import csv
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import IO

from get_file_metadata_helper import get_file_metadata
from pybedtools import BedTool

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths


def mp_process_cnv_data(
    cbio_sample: str,
    cnv_samp_attr: dict,
    info_samp_attr: dict | None,
    seg_samp_attr: dict | None,
    ref_bed: BedTool,
    config_data: dict,
) -> tuple[dict[str, int], dict[str, int], int, list[str], str]:
    """Multi process CNV data by project."""
    ploidy: int = 2
    if cnv_samp_attr["manifest_ftype"] == "ctrlfreec_pval":
        print(f"Processing {cbio_sample} ControlFreeC pval file", file=sys.stderr)

        if info_samp_attr:
            info_fname: str = f"{info_samp_attr['manifest_ftype']}/{info_samp_attr['fname']}"
            ploidy: int = get_ctrlfreec_ploidy(info_fname)
        else:
            print(
                f"WARNING: No info file for {cbio_sample}, using default ploidy of 2",
                file=sys.stderr,
            )
        if seg_samp_attr:
            seg_fname: str = f"{seg_samp_attr['manifest_ftype']}/{seg_samp_attr['fname']}"
        else:
            print(
                f"ERROR: No seg file for {cbio_sample}, cannot process CNV data",
                file=sys.stderr,
            )
            sys.exit(1)

        out_seg_list = read_and_process_ctrlfreec_seg(
            ctrlfreec_seg_fname=seg_fname, sample_id=cbio_sample
        )
        pval_fname: str = f"{cnv_samp_attr['manifest_ftype']}/{cnv_samp_attr['fname']}"
        cnv_bed_obj: BedTool = read_and_process_ctrlfreec_pval(pval_fname=pval_fname)
        raw_cnv_dict, gistic_cnv_dict = get_gene_cnv_dict(
            cnv_bed_obj=cnv_bed_obj,
            ref_bed=ref_bed,
            ploidy=ploidy,
            high_gain=config_data["cnv_high_gain"],
        )
    else:
        print(f"Processing {cbio_sample} GATK CNV file", file=sys.stderr)
        gatk_cnv_fname: str = f"{cnv_samp_attr['manifest_ftype']}/{cnv_samp_attr['fname']}"
        cnv_bed_obj, out_seg_list = read_and_process_gatk_cnv(
            gatk_seg_fname=gatk_cnv_fname, sample_id=cbio_sample
        )
        raw_cnv_dict, gistic_cnv_dict = get_gene_cnv_dict(
            cnv_bed_obj=cnv_bed_obj, ref_bed=ref_bed, high_gain=config_data["cnv_high_gain"]
        )

    return raw_cnv_dict, gistic_cnv_dict, ploidy, out_seg_list, cbio_sample


def output_wide_format_table(
    cnv_dict: dict[str, dict], out_filename: str, ploidy_dict: dict[str, int], cnv_type: str
) -> None:
    """Output wide format table of CNV data.

    Args:
        cnv_dict: Dictionary with gene names as keys and CNV values as values
        out_filename: Output filename for wide format table
        ploidy_dict: Dictionary with sample names as keys and ploidy values as values
        cnv_type: Type of CNV data ("raw" or "gistic")

    """
    # Extract all unique genes from the dictionary
    sample_names = list(cnv_dict.keys())
    all_genes = set()
    for sample in sample_names:
        all_genes.update(cnv_dict[sample].keys())
    all_genes = sorted(all_genes)

    # Open the output file for writing
    with open(out_filename, "w") as out_file:
        # Write the header row with sample names
        print("Hugo_Symbol", *sample_names, file=out_file, sep="\t")
        # Write the data rows for each gene
        for gene in all_genes:
            # Create a list to hold the CNV values for each sample
            cnv_values = []
            if cnv_type == "raw":
                for sample in sample_names:
                    # Get the CNV value for the gene in the current sample, or ploidy if not present
                    cnv_value = cnv_dict[sample].get(gene, ploidy_dict[sample])
                    cnv_values.append(cnv_value)
            else:
                for sample in sample_names:
                    # Get the CNV value for the gene in the current sample, or 0 if not present
                    cnv_value = cnv_dict[sample].get(gene, 0)
                    cnv_values.append(cnv_value)
            # Write the gene name and CNV values to the output file
            print(*[gene, *cnv_values], file=out_file, sep="\t")


def get_gene_cnv_dict(
    cnv_bed_obj: BedTool, ref_bed: BedTool, ploidy: int = 2, high_gain: int = 4
) -> tuple[dict[str, int], dict[str, int]]:
    """Annotate CNV data with gene names and convert to GISTIC-style values.

    Args:
        cnv_bed_obj: BedTool object with CNV data
        ref_bed: BedTool object with reference gene data
        ploidy: Ploidy value for the sample
        high_gain: High gain threshold for GISTIC-style values

    Returns:
        raw_cnv_dict: Dictionary with gene names as keys and CNV values as values

    """
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
        if (
            gene not in raw_cnv_dict
            or (size > raw_cnv_dict[gene][0] and cn != ploidy)
            or (
                size == raw_cnv_dict[gene][0]
                and abs(ploidy - cn) > abs(ploidy - raw_cnv_dict[gene][1])
            )
        ):
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


def read_and_process_ctrlfreec_seg(
    ctrlfreec_seg_fname: str, sample_id: str
) -> list[str]:
    """Process ControlFreeC seg file and write to output.

    Replaces sample ID with cBio sample ID and strips leading chr from chromosome name, then writes to output file.
    Args:
        ctrlfreec_seg_fname: Filename of ControlFreeC seg file
        out_seg_io: Output cbio seg file handle
        sample_id: Sample ID to use

    """
    with open(ctrlfreec_seg_fname) as f:
        cnv_reader: csv.reader._reader = csv.reader(f, delimiter="\t")
        _header = next(cnv_reader)
        out_seg_list: list[str] = []
        for entry in cnv_reader:
            entry[0] = sample_id
            entry[1] = entry[1][3:]  # remove leading chr
            out_seg_list.append("\t".join(entry))
        return out_seg_list


def get_ctrlfreec_ploidy(cfree_info_fname: str) -> int:
    """Process ControFreeC Info file for Ploidy information.

    Used to set a default value for CN if sample does not have a one for a specific gene
    Args:
        cfree_info_fname: Filename of ControFreeC info file
    Returns: int of output ploidy
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


def read_and_process_ctrlfreec_pval(pval_fname: str) -> BedTool:
    """Process ControlFreeC pval file and return BedTool object.

    Filter CNV data by pvalues, strip leading chr, and convert to BedTool object

    Args:
        pval_fname: Filename of ControlFreeC pval file
    Returns: BedTool object with filtered CNV data

    """
    PVALUE: float = 0.05
    with open(pval_fname) as f:
        cnv_reader: csv.reader._reader = csv.reader(f, delimiter="\t")
        header = next(cnv_reader)
        wilcox_idx: int = header.index("WilcoxonRankSumTestPvalue")
        ks_idx: int = header.index("KolmogorovSmirnovPvalue")
        # create bed str representation of CNV data, filtering by pvalue
        cnv_filtered_as_bed: str = "\n".join(
            [
                "\t".join(bed[0:4])
                for bed in cnv_reader
                if (bed[wilcox_idx] != "NA" and bed[ks_idx] != "NA")
                and (float(bed[wilcox_idx]) <= PVALUE and float(bed[ks_idx]) <= PVALUE)
            ]
        )
        cnv_bed_obj: BedTool = BedTool(cnv_filtered_as_bed, from_string=True)
        return cnv_bed_obj


def read_and_process_gatk_cnv(
    gatk_seg_fname: str, sample_id: str
) -> tuple[BedTool, list[str]]:
    """Process GATK CNV seg file and return BedTool object.

    Convert MEAN_LOG2_COPY_RATIO to CN and strip leading chr from chromosome name
    Also write to seg file as it is processed since ibnput is already a seg file

    Args:
        gatk_seg_fname: Filename of GATK CNV seg file
        sample_id: Sample ID to use

    Returns:
        BedTool object with CNV data
        List of seg file entries to write to output file

    """
    with open(gatk_seg_fname) as f:
        cnv_reader: csv.reader._reader = csv.reader(f, delimiter="\t")
        # Need to skip ridiculous interval list-style headers
        cnv_list = [entry for entry in cnv_reader if not entry[0].startswith("@")]
        # Drop chrM, convert MEAN_LOG2_COPY_RATIO to CN and rm leading chr
        # Also output to merged seg file while in the process
        cnv_as_bed = ""
        out_seg_list: list[str] = []
        for entry in cnv_list[1:]:
            if entry[0] != "chrM":
                entry[0] = entry[0][3:]
                out_seg_list.append("\t".join([sample_id] + entry[0:4]))
                # ratio is typically log2(cn/ploidy)
                entry[-2] = round(pow(2, float(entry[-2])) * 2)
                cnv_as_bed += f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[-2]}\n"
        cnv_bed_obj: BedTool = BedTool(cnv_as_bed, from_string=True)
        return cnv_bed_obj, out_seg_list


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
    # cnv meta should have ALL IDs to be processed
    for project, samp_data in cnv_meta.items():
        # initialize out_seg_io here since it can be appended on-the-fly
        out_seg_io: IO = open(f"{out_dir}/{project}.merged_seg.txt", "w")
        print("ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean", file=out_seg_io)
        # store raw and GISTIC-transformed CNV datato be transposed into wide format
        raw_cnv_dict: dict[str, dict] = {}
        gistic_cnv_dict: dict[str, dict] = {}
        # Keep track of ploidy for ControlFreeC to use as default value
        ploidy_dict: dict[str, int] = {}
        samp_list: list[str] = list(samp_data.keys())
        with ProcessPoolExecutor() as executor:
            tasks = [
                executor.submit(
                    mp_process_cnv_data,
                    cbio_sample,
                    cnv_meta[project][cbio_sample],
                    (
                        info_meta[project][cbio_sample] # ControlFreeC specific
                        if project in info_meta and cbio_sample in info_meta[project]
                        else None
                    ),
                    (
                        seg_meta[project][cbio_sample] # ControlFreeC specific
                        if project in seg_meta and cbio_sample in seg_meta[project]
                        else None
                    ),
                    ref_bed,
                    config_data,
                )
                for cbio_sample in samp_list
            ]
            for task in as_completed(tasks):
                ss_raw_cnv_dict, ss_gistic_cnv_dict, ploidy, out_seg_list, cbio_sample = task.result()
                raw_cnv_dict[cbio_sample] = ss_raw_cnv_dict
                gistic_cnv_dict[cbio_sample] = ss_gistic_cnv_dict
                ploidy_dict[cbio_sample] = ploidy
                print(*out_seg_list, sep="\n", file=out_seg_io)
        out_seg_io.close()
        # Print raw and GISTIC results to wide format
        print(f"Writing {project} CNV data to wide format", file=sys.stderr)
        output_wide_format_table(
            cnv_dict=raw_cnv_dict,
            out_filename=f"{out_dir}/{project}.predicted_cnv.txt",
            ploidy_dict=ploidy_dict,
            cnv_type="raw",
        )
        output_wide_format_table(
            cnv_dict=gistic_cnv_dict,
            out_filename=f"{out_dir}/{project}.discrete_cnvs_cnv.txt",
            ploidy_dict=ploidy_dict,
            cnv_type="gistic",
        )


if __name__ == "__main__":
    main()
