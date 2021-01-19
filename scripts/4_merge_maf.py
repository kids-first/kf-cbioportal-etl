#!/usr/bin/env python3
import sys
import os
import argparse
import json
import subprocess
from .get_file_metadata_helper import get_file_metadata
import concurrent.futures


def process_maf(maf_fn, new_maf, maf_exc, tum_id, norm_id, tid_idx, nid_idx, v_idx, eid_idx, n_idx):
    cur_maf = open(maf_fn)
    next(cur_maf)
    next(cur_maf)
    for maf in cur_maf:
        data = maf.rstrip('\n').split('\t')
        if data[v_idx] not in maf_exc:
            data[tid_idx] = tum_id
            data[nid_idx] = norm_id
            data.pop(eid_idx)
            data.pop((n_idx-1))
            new_maf.write('\t'.join(data) + '\n')
    cur_maf.close()


def process_tbl(maf_dir, cbio_dx, file_meta_dict, tid_idx, nid_idx, v_idx, eid_idx, n_idx, print_head):
    try:
        x = 0
        # project/disease name should be name of directory hosting datasheet
        sys.stderr.write('Processing ' + cbio_dx + ' project' + '\n')
        new_maf = open(out_dir + cbio_dx + ".maf", 'w')
        new_maf.write(print_head)
        for cbio_tum_id in file_meta_dict[cbio_dx]:
            cbio_norm_id = file_meta_dict[cbio_dx][cbio_tum_id]['cbio_norm_id']
            fname = file_meta_dict[cbio_dx][cbio_tum_id]['fname']
            sys.stderr.write('Found relevant maf to process for ' + ' ' + cbio_tum_id + ' ' + cbio_norm_id + ' '
            + file_meta_dict[cbio_dx][cbio_tum_id]['kf_tum_id'] + ' ' + file_meta_dict[cbio_dx][cbio_tum_id]['kf_norm_id'] + ' ' + fname + '\n')
            sys.stderr.flush()
            process_maf(maf_dir + fname, new_maf, maf_exc, cbio_tum_id, cbio_norm_id, tid_idx, nid_idx, v_idx, eid_idx, n_idx)
            x += 1
        sys.stderr.write('Completed processing ' + str(x) + ' entries in ' + cbio_dx + '\n')
        new_maf.close()
    except Exception as e:
        print(e)
        sys.exit()

def prepareMAF(config_data):
    # get maf file ext
    maf_dir = '/Users/kalletlak/Documents/temorary/data'
    if maf_dir[-1] != '/':
        maf_dir += '/'
    file_meta_dict = get_file_metadata(args.table, 'maf')

    cur_head = 'Hugo_Symbol	Entrez_Gene_Id	Center	NCBI_Build	Chromosome	Start_Position	End_Position	Strand	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	dbSNP_RS	dbSNP_Val_Status	Tumor_Sample_Barcode	Matched_Norm_Sample_Barcode	Match_Norm_Seq_Allele1	Match_Norm_Seq_Allele2	Tumor_Validation_Allele1	Tumor_Validation_Allele2	Match_Norm_Validation_Allele1	Match_Norm_Validation_Allele2	Verification_Status	Validation_Status	Mutation_Status	Sequencing_Phase	Sequence_Source	Validation_Method	Score	BAM_File	Sequencer	Tumor_Sample_UUID	Matched_Norm_Sample_UUID	HGVSc	HGVSp	HGVSp_Short	Transcript_ID	Exon_Number	t_depth	t_ref_count	t_alt_count	n_depth	n_ref_count	n_alt_count	all_effects	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	ALLELE_NUM	DISTANCE	STRAND_VEP	SYMBOL	SYMBOL_SOURCE	HGNC_ID	BIOTYPE	CANONICAL	CCDS	ENSP	SWISSPROT	TREMBL	UNIPARC	RefSeq	SIFT	PolyPhen	EXON	INTRON	DOMAINS	AF	AFR_AF	AMR_AF	ASN_AF	EAS_AF	EUR_AF	SAS_AF	AA_AF	EA_AF	CLIN_SIG	SOMATIC	PUBMED	MOTIF_NAME	MOTIF_POS	HIGH_INF_POS	MOTIF_SCORE_CHANGE	IMPACT	PICK	VARIANT_CLASS	TSL	HGVS_OFFSET	PHENO	MINIMISED	ExAC_AF	ExAC_AF_AFR	ExAC_AF_AMR	ExAC_AF_EAS	ExAC_AF_FIN	ExAC_AF_NFE	ExAC_AF_OTH	ExAC_AF_SAS	GENE_PHENO	FILTER	flanking_bps	vcf_id	vcf_qual	ExAC_AF_Adj	ExAC_AC_AN_Adj	ExAC_AC_AN	ExAC_AC_AN_AFR	ExAC_AC_AN_AMR	ExAC_AC_AN_EAS	ExAC_AC_AN_FIN	ExAC_AC_AN_NFE	ExAC_AC_AN_OTH	ExAC_AC_AN_SAS	ExAC_FILTER	gnomAD_AF	gnomAD_AFR_AF	gnomAD_AMR_AF	gnomAD_ASJ_AF	gnomAD_EAS_AF	gnomAD_FIN_AF	gnomAD_NFE_AF	gnomAD_OTH_AF	gnomAD_SAS_AF	vcf_pos'
    cur_header = cur_head.rstrip('\n').split('\t')
    eid_idx = cur_header.index('Entrez_Gene_Id')
    n_idx = cur_header.index('NCBI_Build')
    tid_idx = cur_header.index('Tumor_Sample_Barcode')
    nid_idx = cur_header.index('Matched_Norm_Sample_Barcode')
    v_idx = cur_header.index('Variant_Classification')
    cur_header.pop(eid_idx)
    cur_header.pop((n_idx-1))

    print_head = '\t'.join(cur_head) + '\n'
    maf_exc = {"Silent": 0, "Intron": 0, "IGR": 0, "3'UTR": 0, "5'UTR": 0, "3'Flank": 0, "5'Flank": 0, "RNA": 0}
    out_dir = 'merged_mafs/'
    try:
        os.mkdir(out_dir)
    except:
        sys.stderr.write('output dir already exists\n')

    with concurrent.futures.ProcessPoolExecutor(config_data['cpus']) as executor:
        results = {executor.submit(process_tbl, cbio_dx, file_meta_dict, tid_idx, nid_idx, v_idx, eid_idx, n_idx, print_head): cbio_dx for cbio_dx in file_meta_dict}

    sys.stderr.write('Done, check logs\n')
