#!/usr/bin/env python
"""Read in GATK and ControlFreeC CNV data and convert to cBio format."""

import argparse
import csv
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import IO

from pybedtools import BedTool, cleanup

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths


def mp_process_cnv_data(
    cbio_sample: str,
    cnv_samp_attr: dict[str, str],
    info_samp_attr: dict[str, str] | None,
    seg_samp_attr: dict[str, str] | None,
    ref_bed: BedTool,
    config_data: dict,
) -> tuple[dict[str, int], dict[str, int], int, list[str], str]:
    """Multi process CNV data by project."""
    ploidy: int = 2
    try:
        cnv_manifest_ftype: str = cnv_samp_attr["manifest_ftype"]
        if cnv_manifest_ftype == "ctrlfreec_pval" or cnv_manifest_ftype.endswith("_calls"):
            print(f"Processing {cbio_sample} {cnv_manifest_ftype} file", file=sys.stderr)
            if seg_samp_attr:
                seg_fname: str = f"{seg_samp_attr['manifest_ftype']}/{seg_samp_attr['fname']}"
            else:
                print(
                    f"ERROR: No seg file for {cbio_sample}, cannot process CNV data",
                    file=sys.stderr,
                )
                sys.exit(1)
            # ControlFreeC info file is optional, but if it exists, use it to set ploidy
            if info_samp_attr:
                info_fname: str = f"{info_samp_attr['manifest_ftype']}/{info_samp_attr['fname']}"
                ploidy: int = get_ctrlfreec_ploidy(info_fname)
            else:
                print(
                    f"WARNING: No info file for {cbio_sample}, using default ploidy of 2",
                    file=sys.stderr,
                )

            out_seg_list = read_and_process_ctrlfreec_seg(
                ctrlfreec_seg_fname=seg_fname, sample_id=cbio_sample
            )
            cnv_fname: str = f"{cnv_manifest_ftype}/{cnv_samp_attr['fname']}"
            # only a slight difference iv cnvkit and pval inputs, the rest should be the same
            if cnv_manifest_ftype == "ctrlfreec_pval":
                cnv_bed_obj: BedTool = read_and_process_ctrlfreec_pval(pval_fname=cnv_fname)
            else:
                cnv_bed_obj: BedTool = read_and_process_cnvkit_theta2_cns(cns_fname=cnv_fname)
            raw_cnv_dict, gistic_cnv_dict = get_gene_cnv_dict(
                cnv_bed_obj=cnv_bed_obj,
                ref_bed=ref_bed,
                ploidy=ploidy,
                high_gain=config_data["cnv_high_gain"],
            )
        else:
            print(f"Processing {cbio_sample} GATK CNV file", file=sys.stderr)
            gatk_cnv_fname: str = f"{cnv_manifest_ftype}/{cnv_samp_attr['fname']}"
            cnv_bed_obj, out_seg_list = read_and_process_gatk_cnv(
                gatk_seg_fname=gatk_cnv_fname, sample_id=cbio_sample
            )
            raw_cnv_dict, gistic_cnv_dict = get_gene_cnv_dict(
                cnv_bed_obj=cnv_bed_obj, ref_bed=ref_bed, high_gain=config_data["cnv_high_gain"]
            )
    except Exception as e:
        print(f"ERROR: {e} processing {cbio_sample}", file=sys.stderr)
        # clear out pybed temp files
        cleanup(remove_all=True)
        sys.exit()
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
    try:
        # Extract all unique genes from the dictionary
        sample_names = cnv_dict.keys()
        all_genes = set()
        for sample in sample_names:
            all_genes.update(cnv_dict[sample].keys())

        # Open the output file for writing
        with open(out_filename, "w") as out_file:
            # Write the header row with sample names
            print("Hugo_Symbol", *sample_names, file=out_file, sep="\t")
            # Write the data rows for each gene
            for gene in all_genes:
                # Create a list to hold the CNV values for each sample
                default_cnv = 0 if cnv_type == "raw" else ploidy_dict[sample]
                # Get the CNV value for the gene in the current sample, or default if not present
                cnv_values = [cnv_dict[sample].get(gene, default_cnv) for sample in sample_names]
                # Write the gene name and CNV values to the output file
                print(gene, *cnv_values, sep="\t", file=out_file)
    except Exception as e:
        print(f"ERROR: {e} writing {out_filename}", file=sys.stderr)
        # clear out pybed temp files
        cleanup(remove_all=True)
        sys.exit(1)


def get_gene_cnv_dict(
    cnv_bed_obj: BedTool, ref_bed: BedTool, ploidy: int = 2, high_gain: int = 4
) -> tuple[dict[str, int], dict[str, int]]:
    """Annotate CNV data with gene names and convert to GISTIC-style values.

    Args:
        cnv_bed_obj: BedTool object with CNV data
        ref_bed: BedTool object with reference gene data
        ploidy: Ploidy value for the sample
        high_gain: Value over ploidy that results in a GISTIC-style high gain designation

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
    # cleanup pybedtools temp files
    annotated.delete_temporary_history(ask=False)
    return raw_cnv_dict, gistic_cnv_dict


def read_and_process_ctrlfreec_seg(ctrlfreec_seg_fname: str, sample_id: str) -> list[str]:
    """Process ControlFreeC seg file and write to output.

    Replaces sample ID with cBio sample ID and strips leading chr from chromosome name, then writes to output file.
    Args:
        ctrlfreec_seg_fname: Filename of ControlFreeC seg file
        sample_id: Sample ID to use

    """
    with open(ctrlfreec_seg_fname) as f:
        cnv_reader: csv.reader._reader = csv.reader(f, delimiter="\t")
        _header = next(cnv_reader)
        out_seg_list: list[str] = []
        for entry in cnv_reader:
            entry[0] = sample_id
            entry[1] = entry[1].removeprefix("chr")
            out_seg_list.append("\t".join(entry))
    return out_seg_list


def get_ctrlfreec_ploidy(cfree_info_fname: str) -> int:
    """Process ControFreeC Info file for Ploidy information.

    Used to set a default value for CN if sample does not have a one for a specific gene
    Args:
        cfree_info_fname: Filename of ControFreeC info file
    Returns: int of output ploidy
    """
    ploidy = 2
    with open(cfree_info_fname) as info:
        for data in info:
            datum: str = data.rstrip("\n")
            if datum.startswith("Output_Ploidy"):
                try:
                    (_key, value) = datum.split("\t")
                    ploidy = int(value)
                except Exception as e:
                    print(f"WARN: {e} entry could not be split", file=sys.stderr)
                break
    return ploidy


def read_and_process_ctrlfreec_pval(pval_fname: str) -> BedTool:
    """Process ControlFreeC pval file and return BedTool object.

    Filter CNV data by pvalues, convert to BedTool object

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


def read_and_process_cnvkit_theta2_cns(cns_fname: str) -> BedTool:
    """Process CNS file from CNVkit/theta2 file and return BedTool object.

    Strip leading chr before converting to BedTool object

    Args:
        cns_fname: Filename of CNS file
    Returns: BedTool object with filtered CNV data

    """
    with open(cns_fname) as f:
        cnv_reader = csv.reader(f, delimiter="\t")
        _header = next(cnv_reader)
        # create bed str representation of CNV data, stripping leading chr from chrom
        cnv_filtered_as_bed: str = "\n".join(
            [
                "\t".join([bed[0].removeprefix("chr"), bed[1], bed[2], bed[6]])
                for bed in cnv_reader
            ]
        )
        cnv_bed_obj: BedTool = BedTool(cnv_filtered_as_bed, from_string=True)
    return cnv_bed_obj


def read_and_process_gatk_cnv(gatk_seg_fname: str, sample_id: str) -> tuple[BedTool, list[str]]:
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
        cnv_as_bed = ""
        out_seg_list: list[str] = []
        # Skip interval list-style headers by stopping at CONTIG
        # Skip chrM
        # Strip chr from chromosome names
        # Output to merged seg file while in the process
        # Convert MEAN_LOG2_COPY_RATIO to CN and rm leading chr
        for line in cnv_reader:
            if line[0].startswith("CONTIG"):
                break
        for line in cnv_reader:
            if line[0].startswith("chrM"):
                continue
            line[0] = line[0].removeprefix("chr")
            out_seg_list.append("\t".join([sample_id] + line[0:5]))
            # ratio is typically log2(cn/ploidy)
            line[-2] = round(pow(2, float(line[-2])) * 2)
            cnv_as_bed += f"{line[0]}\t{line[1]}\t{line[2]}\t{line[-2]}\n"
        cnv_bed_obj: BedTool = BedTool(cnv_as_bed, from_string=True)
    return cnv_bed_obj, out_seg_list


def pioritize_cnvs(
    cbio_tbl, config_data: dict[str, list[str]]
) -> dict[str, dict[str, dict[str, dict[str, str]]]]:
    """Prioritize CNVs based on config data.

    Args:
        cbio_tbl: CSV reader object with cBio etl TSV
        config_data: CNVs defined by cnv_priority

    Returns:
        Dict of cnv types with prioritized file names and manifest file types.

    """
    # lists and dicts of CNV types
    cnv_ftypes: dict[str, dict[str, dict[str, dict[str, str]]]] = {"cnv": {}, "seg": {}, "info": {}}

    header = next(cbio_tbl)
    # get key fields from header
    cbio_project_idx = header.index("cbio_project")
    cbio_sample_idx = header.index("cbio_sample_name")
    etl_ftype_idx = header.index("etl_file_type")
    etl_exp_idx = header.index("etl_experiment_strategy")
    manifest_ftype_idx = header.index("file_type")
    fname_idx = header.index("file_name")
    for entry in cbio_tbl:
        project = entry[cbio_project_idx]
        cbio_sample = entry[cbio_sample_idx]
        etl_ftype = entry[etl_ftype_idx]
        etl_exp = entry[etl_exp_idx]
        manifest_ftype = entry[manifest_ftype_idx]
        fname = entry[fname_idx]

        if etl_ftype not in cnv_ftypes:
            continue
        if etl_ftype == "info":
            if project not in cnv_ftypes[etl_ftype]:
                cnv_ftypes[etl_ftype][project] = {}
            if cbio_sample not in cnv_ftypes[etl_ftype][project]:
                cnv_ftypes[etl_ftype][project][cbio_sample] = {}
            cnv_ftypes[etl_ftype][project][cbio_sample]["fname"] = fname
        else:
            if project not in cnv_ftypes[etl_ftype]:
                cnv_ftypes[etl_ftype][project] = {}
            if cbio_sample not in cnv_ftypes[etl_ftype][project]:
                cnv_ftypes[etl_ftype][project][cbio_sample] = {}
                cnv_ftypes[etl_ftype][project][cbio_sample]["fname"] = fname
                cnv_ftypes[etl_ftype][project][cbio_sample]["manifest_ftype"] = manifest_ftype
            # If the sample already exists, check if the file type is higher priority
            elif config_data[etl_exp].index(manifest_ftype) < config_data[etl_exp].index(
                cnv_ftypes[etl_ftype][project][cbio_sample]["manifest_ftype"]
            ):
                # If the new file type is higher priority, update the sample entry
                cnv_ftypes[etl_ftype][project][cbio_sample]["fname"] = fname
                cnv_ftypes[etl_ftype][project][cbio_sample]["manifest_ftype"] = manifest_ftype

    return cnv_ftypes


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

    # subset cnv and seg data using priority from config file
    with open(args.table) as f:
        cbio_reader: csv.reader._reader = csv.reader(f, delimiter="\t")
        prioritized_cnv_meta = pioritize_cnvs(cbio_reader, config_data["cnv_priority"])

    ref_bed = BedTool(config_data["bed_genes"])
    out_dir: str = "merged_cnvs"
    os.makedirs(out_dir, exist_ok=True)
    # cnv meta should have ALL IDs to be processed
    for project, samp_data in prioritized_cnv_meta["cnv"].items():
        # initialize out_seg_io here since it can be appended on-the-fly
        out_seg_io: IO = open(f"{out_dir}/{project}.merged_seg.txt", "w")
        print("ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean", file=out_seg_io)
        # store raw and GISTIC-transformed CNV datato be transposed into wide format
        raw_cnv_dict: dict[str, dict] = {}
        gistic_cnv_dict: dict[str, dict] = {}
        # Keep track of ploidy for non-GATK to use as default value
        ploidy_dict: dict[str, int] = {}
        samp_list: list[str] = list(samp_data.keys())
        with ProcessPoolExecutor() as executor:
            tasks = [
                executor.submit(
                    mp_process_cnv_data,
                    cbio_sample,
                    prioritized_cnv_meta["cnv"][project][cbio_sample],
                    (
                        prioritized_cnv_meta["info"][project][cbio_sample]  # Non-GATK specific
                        if project in prioritized_cnv_meta["info"]
                        and cbio_sample in prioritized_cnv_meta["info"][project]
                        else None
                    ),
                    (
                        prioritized_cnv_meta["seg"][project][cbio_sample]  # Non-GATK specific
                        if project in prioritized_cnv_meta["seg"]
                        and cbio_sample in prioritized_cnv_meta["seg"][project]
                        else None
                    ),
                    ref_bed,
                    config_data,
                )
                for cbio_sample in samp_list
            ]
            try:
                for task in as_completed(tasks):
                    ss_raw_cnv_dict, ss_gistic_cnv_dict, ploidy, out_seg_list, cbio_sample = (
                        task.result()
                    )
                    raw_cnv_dict[cbio_sample] = ss_raw_cnv_dict
                    gistic_cnv_dict[cbio_sample] = ss_gistic_cnv_dict
                    ploidy_dict[cbio_sample] = ploidy
                    print(*out_seg_list, sep="\n", file=out_seg_io)
            except KeyboardInterrupt:
                print("Keyboard interrupt, exiting", file=sys.stderr)
                executor.shutdown(wait=False, cancel_futures=True)
                sys.exit(1)
            except TimeoutError:
                print("Timeout occurred during shutdown.", file=sys.stderr)
            finally:
                cleanup(remove_all=True)
                out_seg_io.close()
        # Print raw and GISTIC results to wide format
        print(f"Writing {project} raw CNV data to wide format", file=sys.stderr)
        output_wide_format_table(
            cnv_dict=raw_cnv_dict,
            out_filename=f"{out_dir}/{project}.predicted_cnv.txt",
            ploidy_dict=ploidy_dict,
            cnv_type="raw",
        )
        print(f"Writing {project} GISTIC CNV data to wide format", file=sys.stderr)
        output_wide_format_table(
            cnv_dict=gistic_cnv_dict,
            out_filename=f"{out_dir}/{project}.discrete_cnvs.txt",
            ploidy_dict=ploidy_dict,
            cnv_type="gistic",
        )
        print(f"Finished processing {project} CNV data", file=sys.stderr)


if __name__ == "__main__":
    main()
