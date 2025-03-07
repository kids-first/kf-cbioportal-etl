#!/usr/bin/env python3
import argparse
import concurrent.futures
import json
import os
import subprocess
import sys

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths


def process_cnv(cpath: str) -> None:
    """Process CNVs from a ControlFreeC pvalue output file to gene level.

    Calls bedtools to assign a gene name, and filters out results with < 0.05 P values
    Args:
    cpath: Path of ControlFreeC CNV file
    config_data: dict inherited frm global scope with tool anf referene file locations
    """
    PVALUE = 0.05
    try:
        if cpath != "":
            bedtools = config_data["bedtools"]
            cnv_cp_only = config_data["cp_only_script"]
            bed_file = config_data["bed_genes"]
            hugo_tsv = config_data["hugo_tsv"]
            print("Processing {cpath}", file=sys.stderr)
            root = os.path.basename(cpath)
            temp_filtered_fname = root + ".CNVs_0.05_filtered.bed"
            with open(out_dir + temp_filtered_fname, "w") as temp_filtered, open(cpath) as in_cnv:
                head = next(in_cnv)
                header = head.rstrip("\n").split("\t")
                wilcox_idx = header.index("WilcoxonRankSumTestPvalue")
                ks_idx = header.index("KolmogorovSmirnovPvalue")
                for cnv in in_cnv:
                    ct_dict["total_cnvs"] += 1
                    cnv_data = cnv.rstrip("\n").split("\t")
                    if cnv_data[wilcox_idx] == "NA" or cnv_data[ks_idx] == "NA":
                        ct_dict["NA"] += 1
                    elif float(cnv_data[wilcox_idx]) < PVALUE and float(cnv_data[ks_idx]) < PVALUE:
                        print("\t".join(cnv_data[0:5]), file=temp_filtered)
                    else:
                        ct_dict["p_value_filter"] += 1
            temp_genes_fname = root + ".CNVs.Genes"
            to_genes_cmd = f"{bedtools} intersect -a {out_dir}{temp_filtered_fname} -b {bed_file} -wb > {out_dir}{temp_genes_fname}"
            subprocess.call(to_genes_cmd, shell=True)
            out_gene_cnv_only = f"perl {cnv_cp_only} {hugo_tsv} {out_dir}{temp_genes_fname} > {out_dir}{temp_genes_fname}"
            check = subprocess.call(out_gene_cnv_only, shell=True)
            if check:
                print("Error using to genes script. Check log", file=sys.stderr)
                sys.exit(1)
            rm_tmp = f"rm {out_dir}{temp_genes_fname} {out_dir}{temp_filtered_fname}"
            subprocess.call(rm_tmp, shell=True)
    except Exception as e:
        print(f"{e}", file=sys.stderr)
        sys.exit(1)


parser = argparse.ArgumentParser(
    description="Convert ControlFreeC cnv genome coords to gene level and filter based on p value < 0.05"
)
parser.add_argument(
    "-d",
    "--cnv-dir",
    action="store",
    dest="cnv_dir",
    help="cnv as genome coords file directory",
)
parser.add_argument(
    "-j",
    "--config",
    action="store",
    dest="config_file",
    help="json config file with data types and data locations",
)

args = parser.parse_args()
TOOL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

with open(args.config_file) as f:
    config_data = json.load(f)
config_data = resolve_config_paths(config_data, TOOL_DIR)


suffix = config_data["dna_ext_list"]["copy_number"]
print("Getting cnv list", file=sys.stderr)
flist = [
    os.path.join(args.cnv_dir, f)
    for f in os.listdir(args.cnv_dir)
    if os.path.isfile(os.path.join(args.cnv_dir, f))
]
out_dir = "converted_cnvs/"
try:
    os.mkdir(out_dir)
except:
    print("output dir already exists", file=sys.stderr)

ct_dict = {"total_cnvs": 0, "p_value_filter": 0, "NA": 0}
with concurrent.futures.ThreadPoolExecutor(config_data["cpus"]) as executor:
    results = {executor.submit(process_cnv, cpath): cpath for cpath in flist}
print(f"{ct_dict['total_cnvs']} total cnv calls processed", file=sys.stderr)
print(
    f"{ct_dict['p_value_filter']} calls dropped for not meeting p value cutoff of 0.05, {ct_dict['NA']} calls dropped for having an NA p value",
    file=sys.stderr,
)
