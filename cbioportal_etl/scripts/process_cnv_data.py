#!/usr/bin/env python
"""Read in GATK and ControlFreeC CNV data and convert to cBio format."""

import sys
import pdb
import csv
from pybedtools import BedTool


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
        with open("test_raw_cnv.txt", "w") as raw_out, open("test_gistic_cnv.txt", "w") as gistic_out:
            for gene in raw_cnv_dict:
                print(f"{gene}\t{raw_cnv_dict[gene]}", file=raw_out)
                print(f"{gene}\t{gistic_cnv_dict[gene]}", file=gistic_out)


def read_and_process_gatk_cnv(gatk_seg_fname: str, ref_bed: BedTool):
    with open(gatk_seg_fname) as f:
        cnv_reader: csv.reader._reader = csv.reader(f, delimiter="\t")
        # Need to skip ridiculous interval list-style headers
        cnv_list = [entry for entry in cnv_reader if not entry[0].startswith("@")]
        # Drop chrM, convert MEAN_LOG2_COPY_RATIO to CN and rm leading chr
        cnv_as_bed = ""
        for entry in cnv_list[1:]:
            if entry[0] != "chrM":
                entry[0] = entry[0][3:]
                # ratio is typically log2(cn/ploidy)
                delme = entry[-2]
                entry[-2] = pow(2, float(entry[-2])) * 2
                # delme sanity check
                if round(entry[-2]) == 2 and entry[5] == "+":
                    print (*[*entry, delme], sep="\t")
                cnv_as_bed += f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[-2]}\n"
        cnv_bed_obj: BedTool = BedTool(cnv_as_bed, from_string=True)
        raw_cnv_dict, gistic_cnv_dict = get_gene_cnv_dict(cnv_bed_obj, ref_bed)
        with open("gatk_raw_cnv.txt", "w") as raw_out, open("gatk_gistic_cnv.txt", "w") as gistic_out:
            for gene in raw_cnv_dict:
                print(f"{gene}\t{raw_cnv_dict[gene]}", file=raw_out)
                print(f"{gene}\t{gistic_cnv_dict[gene]}", file=gistic_out)

def main():
    ref_bed = BedTool(sys.argv[2])
    # read_and_process_ctrlfreec_pval(sys.argv[1], ref_bed)
    read_and_process_gatk_cnv(sys.argv[1], ref_bed)
if __name__ == "__main__":
    main()
