#!/usr/bin/env python
"""Read in exported germline table, convert to MAF, add to existing somatic MAF."""

import _csv
import argparse
import csv
import json
import os
import sys

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths


def append_variant_type_field(ref_allele, alt_allele):
    """Resolve the VCF variant type.

    Stolen from https://github.com/genome-nexus/vcf2maf-lite/blob/main/vcf2maf_lite/vcf2maf_lite.py#L735
    """
    variant_type = ""

    # first check if indel
    ref_len = len(ref_allele)
    alt_len = len(alt_allele)
    if ref_allele == "-" or ref_len < alt_len:
        variant_type = "INS"
    elif alt_allele == "-" or alt_len < ref_len:
        variant_type = "DEL"
    # check whether variant type is type of polymorphism
    elif ref_len == alt_len:
        if ref_len == 1:
            variant_type = "SNP"
        elif ref_len == 2:
            variant_type = "DNP"
        elif ref_len == 3:
            variant_type = "TNP"
        else:
            variant_type = "ONP"

    # if variant type is still empty then report
    if variant_type == "":
        message = (
            f"Could not salvage variant type from alleles [ ref allele = {ref_allele} , tumor allele = {alt_allele} ]"
        )
        print(message, file=sys.stderr)
    inframe = abs(ref_len - alt_len) % 3 == 0

    return variant_type, inframe


def append_variant_class(effect: str, var_type: str, inframe: int) -> str:
    """Convert VEP consequence to cBioPortal variant classification.

    Args:
        effect (str): VEP consequence term
        var_type (str): Variant type (DEL, INS, SNP)
        inframe (int): 1 if inframe, 0 if not

    Returns:
        str: cBioPortal variant classification

    """
    if not effect:
        return "Targeted_Region"

    if effect == "frameshift_variant" or (
        effect == "protein_altering_variant" and not inframe
    ):
        return {
            "DEL": "Frame_Shift_Del",
            "INS": "Frame_Shift_Ins",
        }.get(var_type, "Targeted_Region")

    if effect == "protein_altering_variant" and inframe:
        return {
            "DEL": "In_Frame_Del",
            "INS": "In_Frame_Ins",
        }.get(var_type, "Targeted_Region")

    if effect.endswith("inframe_insertion"):
        return "In_Frame_Ins"

    if effect.endswith("inframe_deletion"):
        return "In_Frame_Del"

    for classification, effects in EFFECT_MATCHES.items():
        if effect in effects:
            return classification

    return "Targeted_Region"


def append_amino_acids_field(hgvs_p: str) -> str:
    """Convert 3-letter amino acid codes to 1-letter codes."""
    hgvs_p_short = hgvs_p
    for find, replace in AA_3_TO_1.items():
        hgvs_p_short = hgvs_p_short.replace(find, replace)

    return hgvs_p_short


def main() -> None:
    """Parse arguments and process germline MAF."""
    parser = argparse.ArgumentParser(
        description="Clean up and append germline MAF to somatic MAF for cBioPortal import"
    )
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        help="Table with cbio project, kf bs ids, cbio IDs, and file names",
    )

    parser.add_argument(
        "-j",
        "--config",
        action="store",
        dest="config_file",
        help="json config file with data types and data locations",
    )

    parser.add_argument(
        "-e",
        "--exported",
        action="store",
        help="exported germline variants TSV",
    )

    args = parser.parse_args()
    TOOL_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    germ_data_dict: dict[str, list] = {}
    with open(args.exported) as f:
        germline_reader: csv.reader._reader = csv.reader(f, delimiter="\t")
        header = next(germline_reader)
        header.extend(NEW_FIELDS)
        print("\t".join(header))
        # get indices of fields used to inform fields to be calculated
        sample_id_idx = header.index("Matched_Norm_Sample_Barcode")
        ref_idx = header.index("Reference_Allele")
        alt_idx = header.index("Match_Norm_Seq_Allele1")
        hgvs_p_idx = header.index("HGVSp")
        csq_idx = header.index("Consequence")
        for data in germline_reader:
            # add new fields to data
            var_type, inframe = append_variant_type_field(data[ref_idx], data[alt_idx])
            data.append(append_amino_acids_field(data[hgvs_p_idx]))
            data.append(var_type)
            data.append(append_variant_class(data[csq_idx], var_type, inframe))
            # add to dict by sample id
            sample_id = data[sample_id_idx]
            if sample_id not in germ_data_dict:
                germ_data_dict[sample_id] = []
            germ_data_dict[sample_id].append(data)
            # debug, print out data
            print("\t".join(data))

    # with open(args.config_file) as f:
    #     config_data = json.load(f)
    # config_data: dict = resolve_config_paths(config_data, TOOL_DIR)


if __name__ == "__main__":
    # set global vars here, then call main
    AA_3_TO_1 = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Asx": "B",
        "Cys": "C",
        "Glu": "E",
        "Gln": "Q",
        "Glx": "Z",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
        "Xxx": "X",
        "Ter":  "*",
    }

    EFFECT_MATCHES = {
        "Splice_Site": {
            "splice_acceptor_variant",
            "splice_donor_variant",
            "transcript_ablation",
            "exon_loss_variant",
        },
        "Nonsense_Mutation": {"stop_gained"},
        "Nonstop_Mutation": {"stop_lost"},
        "Translation_Start_Site": {
            "initiator_codon_variant",
            "start_lost",
        },
        "Missense_Mutation": {
            "missense_variant",
            "coding_sequence_variant",
            "conservative_missense_variant",
            "rare_amino_acid_variant",
        },
        "Intron": {
            "transcript_amplification",
            "intron_variant",
            "INTRAGENIC",
            "intragenic_variant",
        },
        "Splice_Region": {"splice_region_variant"},
        "Silent": {
            "incomplete_terminal_codon_variant",
            "synonymous_variant",
            "stop_retained_variant",
            "NMD_transcript_variant",
        },
        "RNA": {
            "mature_miRNA_variant",
            "exon_variant",
            "non_coding_exon_variant",
            "non_coding_transcript_exon_variant",
            "non_coding_transcript_variant",
            "nc_transcript_variant",
        },
        "5'UTR": {
            "5_prime_UTR_variant",
            "5_prime_UTR_premature_start_codon_gain_variant",
        },
        "3'UTR": {"3_prime_UTR_variant"},
        "IGR": {
            "TF_binding_site_variant",
            "regulatory_region_variant",
            "regulatory_region",
            "intergenic_variant",
            "intergenic_region",
        },
        "5'Flank": {"upstream_gene_variant"},
        "3'Flank": {"downstream_gene_variant"},
    }

    NEW_FIELDS = ["HGVSp_Short", "Variant_Type", "Variant_Classification"]
    main()
