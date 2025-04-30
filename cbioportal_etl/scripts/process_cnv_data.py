#!/usr/bin/env python
"""Read in GATK and ControlFreeC CNV data and convert to cBio format."""

import sys
import pdb
import csv
from pybedtools import BedTool



def read_and_process_ctrlfreec_pval(pval_fname: str, ref_bed: BedTool) -> None:
    PVALUE = 0.05
    with open(pval_fname) as f:
        test = f.readlines()
        cnv_list = list(csv.reader(test, delimiter="\t"))
        header = cnv_list[0]
        wilcox_idx = header.index("WilcoxonRankSumTestPvalue")
        ks_idx = header.index("KolmogorovSmirnovPvalue")
        # create bed str representation of CNV data, filtering by pvalue
        bed_str = "\n".join(["\t".join(bed[0:4]) for bed in cnv_list[1:] if float(bed[wilcox_idx]) <= PVALUE and float(bed[ks_idx]) <= PVALUE])
        cnv_bed_obj = BedTool(bed_str, from_string=True)
        # Create staging object that has cnv chrom, cnv start, cnv end, cn
        # ref chr, ref start, ref end, ref gene size overlap
        annotated = cnv_bed_obj.intersect(ref_bed, b=True, wo=True)
        # iterate though annotated bed obj as list of lists to store in dict with uniq gene: CN values
        pdb.set_trace()
        hold = 1

    pass


def read_and_process_gatk_cnv():
    pass

def main():
    ref_bed = BedTool(sys.argv[2])
    read_and_process_ctrlfreec_pval(sys.argv[1], ref_bed)

if __name__ == "__main__":
    main()
