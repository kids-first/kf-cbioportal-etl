#!/usr/bin/env python3
import sys
import os
import argparse
import json
import subprocess
import concurrent.futures


def process_cnv(cpath):
    try:
        if cpath != "":
            bedtools = config_data["bedtools"]
            cnv_cp_only = config_data["cp_only_script"]
            bed_file = config_data["bed_genes"]
            hugo_tsv = config_data["hugo_tsv"]
            print("Processing " + cpath, file=sys.stderr)
            root = os.path.basename(cpath)
            temp_filtered_fname = root + ".CNVs_0.05_filtered.bed"
            with open(out_dir + temp_filtered_fname, "w") as temp_filtered:
                in_cnv = open(cpath)
                head = next(in_cnv)
                header = head.rstrip('\n').split('\t')
                wilcox_idx = header.index('WilcoxonRankSumTestPvalue')
                ks_idx = header.index('KolmogorovSmirnovPvalue')
                for cnv in in_cnv:
                    ct_dict["total_cnvs"] += 1
                    cnv_data = cnv.rstrip("\n").split("\t")
                    if cnv_data[wilcox_idx] == "NA" or cnv_data[ks_idx] == "NA":
                        ct_dict["NA"] += 1
                    elif float(cnv_data[wilcox_idx]) < 0.05 and float(cnv_data[ks_idx]) < 0.05:
                        print("\t".join(cnv_data[0:5]), file=temp_filtered)
                    else:
                        ct_dict["p_value_filter"] += 1
            temp_genes_fname = root + ".CNVs.Genes"
            to_genes_cmd = "{} intersect -a {}{} -b {} -wb > {}{}".format(bedtools, out_dir, temp_filtered_fname, bed_file, out_dir, temp_genes_fname)
            subprocess.call(to_genes_cmd, shell=True)
            out_gene_cnv_only = "perl {} {} {}{} > {}{}.copy_number".format(cnv_cp_only, hugo_tsv, out_dir, temp_genes_fname, out_dir, temp_genes_fname)
            check = subprocess.call(out_gene_cnv_only, shell=True)
            if check:
                print('Error using to genes script. Check log', file=sys.stderr)
                sys.exit(1)
            rm_tmp = "rm {}{} {}{}".format(out_dir, temp_genes_fname, out_dir, temp_filtered_fname)
            subprocess.call(rm_tmp, shell=True)
    except Exception as e:
        print("{}".format(e), file=sys.stderr)
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
    help="json config file with data types and " "data locations",
)

args = parser.parse_args()
with open(args.config_file) as f:
    config_data = json.load(f)

suffix = config_data["dna_ext_list"]["copy_number"]
print("Getting cnv list", file=sys.stderr)
flist =[os.path.join(args.cnv_dir, f) for f in os.listdir(args.cnv_dir) if os.path.isfile(os.path.join(args.cnv_dir, f))]
out_dir = "converted_cnvs/"
try:
    os.mkdir(out_dir)
except:
    print("output dir already exists", file=sys.stderr)

ct_dict = {"total_cnvs": 0, "p_value_filter": 0, "NA": 0}
with concurrent.futures.ThreadPoolExecutor(config_data["cpus"]) as executor:
    results = { executor.submit(process_cnv, cpath): cpath for cpath in flist }
print("{} total cnv calls processed".format(ct_dict["total_cnvs"]), file=sys.stderr)
print("{} calls dropped for not meeting p value cutoff of 0.05, {} calls dropped for having an NA p value".format(ct_dict["p_value_filter"], ct_dict["NA"]), file=sys.stderr )
