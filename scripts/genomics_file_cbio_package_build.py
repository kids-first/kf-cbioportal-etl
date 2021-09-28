#!/usr/bin/env python3
"""
This is a genomics file + clinical data file to cBio package conversion workflow script.
It is designed to work with standard KF somatic workflow outputs as well as DGD outputs.
Clinical files should have been produced ahead of time, while supporting sample ID file manifest, case_meta config json file, and data json config.
If there are fusion files, they also need s table outlining sequencing center location.
See README for prequisite details.
"""
import sys
import argparse
import json
import subprocess


def check_status(status, data_type, log_name):
    if status:
        sys.stderr.write(
            "Something when wrong while processing the "
            + data_type
            + "! Check the "
            + log_name
            + " for more info\n"
        )
        exit(1)


def process_maf(maf_loc_dict, cbio_id_table, data_config_file, dgd_status):
    """
    Collate and process pbta/kf style mafs. Call append if also adding dgd data
    """
    sys.stderr.write("Processing maf files\n")
    maf_dir = maf_loc_dict["kf"]
    if args.dgd_status == "dgd":
        maf_dir = maf_loc_dict["dgd"]
    maf_header = maf_loc_dict["header"]
    maf_cmd = "{}maf_merge.py -t {} -i {} -m {} -j {} -f {} 2> collate_mafs.log".format(
        script_dir, cbio_id_table, maf_header, maf_dir, data_config_file, dgd_status
    )
    status = 0
    status = subprocess.call(maf_cmd, shell=True)
    if args.dgd_status == "both":
        status += process_append_dgd_maf(maf_loc_dict, cbio_id_table)
    return status


def process_append_dgd_maf(maf_loc_dict, cbio_id_table):
    """
    Append DGD mafs to existing collated kf maf
    """
    sys.stderr.write("Appending DGD to exsting KF maf\n")
    maf_header = maf_loc_dict["header"]
    in_maf_dir = maf_loc_dict["dgd"]
    append_maf = "merged_mafs/" + cbio_study_id + ".maf"
    append_cmd = (
        "{}add_dgd_maf_to_pbta.py -i {} -m {} -t {} >> {} 2> dgd_append_maf.log".format(
            script_dir, maf_header, in_maf_dir, cbio_id_table, append_maf
        )
    )
    status = 0
    status = subprocess.call(append_cmd, shell=True)
    return status


def process_cnv(cnv_loc_dict, data_config_file, cbio_id_table):
    """
    Add gene info to CNV calls, merge into table, and create GISTIC-style output
    """
    sys.stderr.write("Processing CNV calls\n")
    cnv_dir = cnv_loc_dict["pval"]
    gene_annot_cmd = "{}cnv_1_genome2gene.py -d {} -j {} 2> cnv_gene_annot.log".format(
        script_dir, cnv_dir, data_config_file
    )
    status = 0
    status = subprocess.call(gene_annot_cmd, shell=True)
    info_dir = cnv_loc_dict["info"]
    merge_gene_cmd = "{}cnv_2_merge.py -t {} -n converted_cnvs -j {}".format(
        script_dir, cbio_id_table, data_config_file
    )
    if info_dir != "":
        merge_gene_cmd += " -i " + info_dir
    merge_gene_cmd += " 2> merge_cnv_gene.log"
    status += subprocess.call(merge_gene_cmd, shell=True)
    gistic_style_cmd = "{}cnv_3_gistic_style.py -d merged_cnvs -j {} -t {}".format(
        script_dir, data_config_file, cbio_id_table
    )
    if info_dir != "":
        gistic_style_cmd += " -i " + info_dir
    gistic_style_cmd += " 2> merge_cnv_gene.log"
    status += subprocess.call(gistic_style_cmd, shell=True)

    if "seg" in cnv_loc_dict:
        status += process_seg(cnv_loc_dict, cbio_id_table, data_config_file)

    return status


def process_seg(cnv_loc_dict, cbio_id_table, data_config_file):
    """
    Collate and process CNV seg files
    """
    sys.stderr.write("Processing CNV seg calls\n")
    seg_dir = cnv_loc_dict["seg"]
    merge_seg_cmd = "{}cnv_merge_seg.py -t {} -m {} -j {} 2> cnv_merge_seg.log".format(
        script_dir, cbio_id_table, seg_dir, data_config_file
    )
    status = subprocess.call(merge_seg_cmd, shell=True)
    return status


def process_rsem(rsem_dir, cbio_id_table):
    """
    Merge rsem results by FPKM, calculate z-scores
    """
    sys.stderr.write("Processing RNA expression data\n")
    merge_rsem_cmd = "{}rna_merge_rename_expression.py -t {} -r {} 2> rna_merge_rename_expression.log".format(
        cbio_id_table, rsem_dir
    )
    status = subprocess.call(merge_rsem_cmd, shell=True)
    return status


def process_kf_fusion(fusion_dir, cbio_id_table, sq_file):
    """
    Collate and process annoFuse output
    """
    sys.stderr.write("Processing KF fusion calls\n")
    fusion_cmd = (
        script_dir
        + "rna_convert_fusion.py -t "
        + cbio_id_table
        + " -f "
        + fusion_dir
        + " -m annofuse -s "
        + sq_file
        + " 2> rna_convert_fusion.log"
    )
    fusion_cmd = "{}rna_convert_fusion.py -t {} -f {} -m annofuse -s {} 2> rna_convert_fusion.log".format(
        cbio_id_table, fusion_dir, sq_file
    )
    status = subprocess.call(fusion_cmd, shell=True)
    return status


def process_dgd_fusion(cbio_id_table, fusion_dir, dgd_status):
    """
    Collate process DGD fusion output
    """
    dgd_fusion_cmd = (
        script_dir + "add_dgd_fusion.py -t " + cbio_id_table + " -f " + fusion_dir
    )
    dgd_fusion_cmd = "{}add_dgd_fusion.py -t {} -f {}".format(
        script_dir, cbio_id_table, fusion_dir
    )
    if dgd_status == "both":
        sys.stderr.write("Appending DGD fusion calls\n")
        append_fusion = "merged_fusion/" + cbio_study_id + ".fusions.txt"
        dgd_fusion_cmd += " -a >> " + append_fusion
    else:
        sys.stderr.write("Processing DGD fusion calls\n")
    dgd_fusion_cmd += " 2> add_dgd_fusion.log"
    status = subprocess.call(dgd_fusion_cmd, shell=True)
    return status


# This part is commented out for now as there is no good solution yet for files stored in different account - need to be able to run chopaws to get key
# if(args.manifest != None):
#     sys.stderr.write('Download manifest given. Downloading files\n')
#     dl_file_cmd = config_data['script_dir'] + 'get_files_from_manifest.py -m ' + args.manifest + ' -f ' + join(config_data['dl_file_type_list'] + ' -p saml')
#     subprocess.call(dl_file_cmd, shell=True)

parser = argparse.ArgumentParser(
    description="Download files (if needed), collate genomic files, organize load package."
)
parser.add_argument(
    "-t",
    "--table",
    action="store",
    dest="table",
    help="Table with cbio project, kf bs ids, cbio IDs, and file names",
)
parser.add_argument(
    "-m",
    "--manifest",
    action="store",
    dest="manifest",
    help="Download file manifest, if needed",
)
parser.add_argument(
    "-c",
    "--cbio-config",
    action="store",
    dest="cbio_config",
    help="cbio case and meta config file",
)
parser.add_argument(
    "-d",
    "--data-config",
    action="store",
    dest="data_config",
    help="json config file with data types and " "data locations",
)
parser.add_argument(
    "-f",
    "--dgd-status",
    action="store",
    dest="dgd_status",
    help="Flag to determine load will have pbta/kf + dgd(both), kf/pbta only(kf), dgd-only(dgd)",
    default="both",
    const="both",
    nargs="?",
    choices=["both", "kf", "dgd"],
)


args = parser.parse_args()
with open(args.data_config) as f:
    config_data = json.load(f)

with open(args.cbio_config) as f:
    config_meta_case = json.load(f)

cbio_study_id = config_meta_case["study"]["cancer_study_identifier"]
# iterate through config file - file should only have keys related to data to be loaded
script_dir = config_data["script_dir"]
for key in config_meta_case:
    if key.startswith("merged_"):
        data_type = "_".join(key.split("_")[1:])
        if data_type == "mafs":
            exit_status = process_maf(
                config_data["file_loc_defs"]["mafs"],
                args.table,
                args.data_config,
                args.dgd_status,
            )
            check_status(
                exit_status,
                "mafs and/or dgd_append_maf",
                "collate_mafs.log and/or dgd_append_maf.log",
            )
        elif data_type == "cnvs":
            exit_status = process_cnv(config_data["file_loc_defs"]["cnvs"])
            check_status(exit_status, "cnv data", "cnv_* logs")
        elif data_type == "rsem":
            exit_status = process_rsem(config_data["file_loc_defs"]["rsem"], args.table)
            check_status(
                exit_status, "rna expression", "rna_merge_rename_expression.log"
            )
        elif data_type == "fusion":
            # Status both works for...both, only when one is specifically picked should one not be run
            if args.dgd_status != "dgd":
                exit_status = process_kf_fusion(
                    config_data["file_loc_defs"]["fusion"],
                    args.table,
                    config_data["file_loc_defs"]["fusion_sq_file"],
                )
                check_status(exit_status, "kf fusions", "rna_convert_fusion.log")
            if args.dgd_status != "kf":
                exit_status = process_dgd_fusion(
                    args.table,
                    config_data["file_loc_defs"]["dgd_fusion"],
                    args.dgd_status,
                )
                check_status(exit_status, "dgd fusions", "add_dgd_fusion.log")

# Run final package builder script
sys.stderr.write("Creating load packages\n")
pck_cmd = (
    config_data["script_dir"]
    + "organize_upload_packages.py -o processed -c "
    + args.cbio_config
    + " 2> load_package_create.log"
)
exit_status = subprocess.call(pck_cmd, shell=True)
check_status(exit_status, "load package", "load_package_create.log")

# Run cbioportal data validator
sys.stderr.write("Validating load packages\n")
validate = (
    config_data["cbioportal_validator"]
    + " -s processed/"
    + cbio_study_id
    + " -n -v 2> validator.errs > validator.out"
)
exit_status = subprocess.call(validate, shell=True)
if exit_status:
    sys.stderr.write(
        "Validator quit with status "
        + str(exit_status)
        + ". Check validator.errs and validator.out for more info\n"
    )
