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
        "-a",
        "--append",
        action="store",
        help="MAF to append to",
    )

    parser.add_argument(
        "-e",
        "--exported",
        action="store",
        help="exported germline variants TSV",
    )

    args = parser.parse_args()
    TOOL_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Build a map of cBio IDs to assign to matched normal, then append
    matched_ids: dict[str, set] = {}
    # also get header entries to use for output, if not in germline MAF, will output blank
    h_dict: dict[str, int|None] = {}

    with open(args.append) as in_maf:
        _version_line = in_maf.readline()
        maf_reader: csv.reader._reader = csv.reader(in_maf, delimiter="\t")
        maf_header = next(maf_reader)
        maf_tum_id_idx = maf_header.index("Tumor_Sample_Barcode")
        maf_norm_id_idx = maf_header.index("Matched_Norm_Sample_Barcode")
        for row in maf_reader:
            tum_id, norm_id = row[maf_tum_id_idx], row[maf_norm_id_idx]
            if norm_id not in matched_ids:
                matched_ids[norm_id] = set()
            matched_ids[norm_id].add(tum_id)

    with open(args.append, "a") as in_maf, open(args.exported) as f:
        germline_reader: csv.reader._reader = csv.reader(f, delimiter="\t")

        header = next(germline_reader)
        header.extend(NEW_FIELDS)
        for item in maf_header:
            if item in header:
                h_dict[item] = header.index(item)
            else:
                h_dict[item] = None
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
            # now append all matched tumor IDs for this normal to the MAF, fill in fields not in maf_header
            # repeat for each tumor for now
            if data[sample_id_idx] in matched_ids:
                for tum_id in matched_ids[data[sample_id_idx]]:
                    data_reorg = []
                    for field, idx in h_dict.items():
                        if field == "Hugo_Symbol" and data[header.index("Hugo_Symbol")] == "":
                            print(f"Warning: Hugo_Symbol is empty for {data} in germline MAF", file=sys.stderr)
                        elif idx is not None:
                            data_reorg.append(data[idx])
                        elif field == "Tumor_Sample_Barcode":
                            data_reorg.append(tum_id)
                        elif field == "Matched_Norm_Sample_Barcode":
                            data_reorg.append(data[sample_id_idx])
                        elif field == "Tumor_Seq_Allele1":
                            data_reorg.append(data[alt_idx])
                        else:
                            data_reorg.append("")
                    in_maf.write("\t".join(data_reorg) + "\n")

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
