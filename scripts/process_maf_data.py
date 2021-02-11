#!/usr/bin/env python3
import sys
from .meta_file_utils import add_meta_file


def add_maf_file(config_data, study_name, study_dir, study_maf_files):
    cur_head = 'Hugo_Symbol	Entrez_Gene_Id	Center	NCBI_Build	Chromosome	Start_Position	End_Position	Strand	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	dbSNP_RS	dbSNP_Val_Status	Tumor_Sample_Barcode	Matched_Norm_Sample_Barcode	Match_Norm_Seq_Allele1	Match_Norm_Seq_Allele2	Tumor_Validation_Allele1	Tumor_Validation_Allele2	Match_Norm_Validation_Allele1	Match_Norm_Validation_Allele2	Verification_Status	Validation_Status	Mutation_Status	Sequencing_Phase	Sequence_Source	Validation_Method	Score	BAM_File	Sequencer	Tumor_Sample_UUID	Matched_Norm_Sample_UUID	HGVSc	HGVSp	HGVSp_Short	Transcript_ID	Exon_Number	t_depth	t_ref_count	t_alt_count	n_depth	n_ref_count	n_alt_count	all_effects	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	ALLELE_NUM	DISTANCE	STRAND_VEP	SYMBOL	SYMBOL_SOURCE	HGNC_ID	BIOTYPE	CANONICAL	CCDS	ENSP	SWISSPROT	TREMBL	UNIPARC	RefSeq	SIFT	PolyPhen	EXON	INTRON	DOMAINS	AF	AFR_AF	AMR_AF	ASN_AF	EAS_AF	EUR_AF	SAS_AF	AA_AF	EA_AF	CLIN_SIG	SOMATIC	PUBMED	MOTIF_NAME	MOTIF_POS	HIGH_INF_POS	MOTIF_SCORE_CHANGE	IMPACT	PICK	VARIANT_CLASS	TSL	HGVS_OFFSET	PHENO	MINIMISED	ExAC_AF	ExAC_AF_AFR	ExAC_AF_AMR	ExAC_AF_EAS	ExAC_AF_FIN	ExAC_AF_NFE	ExAC_AF_OTH	ExAC_AF_SAS	GENE_PHENO	FILTER	flanking_bps	vcf_id	vcf_qual	ExAC_AF_Adj	ExAC_AC_AN_Adj	ExAC_AC_AN	ExAC_AC_AN_AFR	ExAC_AC_AN_AMR	ExAC_AC_AN_EAS	ExAC_AC_AN_FIN	ExAC_AC_AN_NFE	ExAC_AC_AN_OTH	ExAC_AC_AN_SAS	ExAC_FILTER	gnomAD_AF	gnomAD_AFR_AF	gnomAD_AMR_AF	gnomAD_ASJ_AF	gnomAD_EAS_AF	gnomAD_FIN_AF	gnomAD_NFE_AF	gnomAD_OTH_AF	gnomAD_SAS_AF	vcf_pos'
    cur_header = cur_head.rstrip('\n').split('\t')
    eid_idx = cur_header.index('Entrez_Gene_Id')
    n_idx = cur_header.index('NCBI_Build')
    tid_idx = cur_header.index('Tumor_Sample_Barcode')
    nid_idx = cur_header.index('Matched_Norm_Sample_Barcode')
    v_idx = cur_header.index('Variant_Classification')
    cur_header.pop(eid_idx)
    cur_header.pop((n_idx-1))
    try:
        sys.stderr.write('Processing for maf ' + study_name + ' project' + '\n')
        new_maf = open(study_dir + "/data_mutations_extended.txt", 'w')
        new_maf.write('\t'.join(cur_header) + '\n')
        for resource in study_maf_files:
            t_bs_id = resource.metadata['Kids First Biospecimen ID Tumor']
            n_bs_id = resource.metadata['Kids First Biospecimen ID Normal']
            process_maf_file(config_data['data_files'] + resource.name, new_maf, t_bs_id, n_bs_id, tid_idx, nid_idx, v_idx, eid_idx, n_idx)
        new_maf.close()

        #Add meta file
        add_meta_file(study_name, config_data["metadata"]["mutation"]["meta_attr"], study_dir + '/meta_mutations_extended.txt')

    except Exception as e:
        print(e)
        sys.exit()

def process_maf_file(maf_fn, new_maf, tum_id, norm_id, tid_idx, nid_idx, v_idx, eid_idx, n_idx):
    maf_exc = {"Silent": 0, "Intron": 0, "IGR": 0, "3'UTR": 0, "5'UTR": 0, "3'Flank": 0, "5'Flank": 0, "RNA": 0}
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
