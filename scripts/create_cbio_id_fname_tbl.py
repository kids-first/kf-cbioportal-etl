#!/usr/bin/env python3
import argparse
import sys
import os
import json
import subprocess
from sevenbridges.models.file import File
import sevenbridges as sbg
import numpy as np
import pandas as pd
import concurrent.futures
from scipy import stats

from .process_fusion_data import add_fusion_file


# def process_ds(dsname, cav_dict, out):
#     try:
#         # last item in find will likely be empty after split
#         x = 0
#         parts = dsname.split('/')
#         cbio_proj = parts[-2]
#         # project/disease name should be name of directory hosting datasheet
#         sys.stderr.write('Processing ' + cbio_proj + ' project' + '\n')
#         cur_ds = open(dsname)
#         for i in range(0, 4, 1):
#             next(cur_ds)
#         ds_head = next(cur_ds)
#         ds_header = ds_head.rstrip('\n').split('\t')
#         sa_idx = ds_header.index('SAMPLE_ID')
#         sp_idx = ds_header.index('SPECIMEN_ID')
#         ns_idx = ds_header.index('MATCHED_NORMAL_SPECIMEN_ID')
#         for entry in cur_ds:
#             meta = entry.rstrip('\n').split('\t')
#             check = meta[sp_idx].split(';')
#             for bs_id in check:
#                 # can't tell if RNA or DNA from datasheet, but DNA will be tum + norm bs id, RNA tum + NA
#                 id1 = bs_id + "\t" + meta[ns_idx]
#                 id2 = bs_id + "\tNA"
#                 key = id1
#                 if key in cav_dict:
#                     for ftype in cav_dict[key]:
#                         out.write("\t".join([cbio_proj, key, ftype, meta[sa_idx], meta[ns_idx], cav_dict[key][ftype]]) + "\n")
#                 elif id2 in cav_dict:
#                     key = id2
#                     for ftype in cav_dict[key]:
#                         out.write("\t".join([cbio_proj, key, ftype, meta[sa_idx], 'NA', cav_dict[key][ftype]]) + "\n")
#                 else:
#                     sys.stderr.write('WARNING: ' + id1 + ' nor ' + id2 + ' found in data sheet entry\n' + entry + 'Check inputs!\n')
#     except Exception as e:
#         print(e)
#         sys.exit()

# parser = argparse.ArgumentParser(description='Create temp table with cbio ids, bs ids, and file locations and types for downstream file merging. Should be run one level above datasheets dir.')
# parser.add_argument('-c', '--cavatica', action='store', dest='cav',
#                     help='file with task info from cavatica (see step 1b)')

# args = parser.parse_args()

# cav_dict = {}
# manifest = open(args.cav)
# head = next(manifest)
# for line in manifest:
#     info = line.rstrip('\n').split('\t')
#     bs_ids = info[0] + "\t" + info[1]
#     fnames = info[-2]
#     atype = info[4]
#     cav_dict[bs_ids] = {}
#     # will likely have to add fusion ext soon
#     for fname in fnames.split(','):
#         ftype = 'rsem'
#         if atype == 'DNA':
#             if fname[-3:] == 'maf':
#                 ftype = 'maf'
#             else:
#                 ftype =  "cnv"
#         elif fname[-3:] == 'tsv':
#             ftype = 'fusion'
#         cav_dict[bs_ids][ftype] = fname

# flist = subprocess.check_output('find /Users/kalletlak/temorary/datasets -name data_clinical_sample.txt', shell=True)

# ds_list = flist.decode().split('\n')
# if ds_list[-1] == '':
#     ds_list.pop()
# out = open('cbio_id_fname_table.txt', 'w')
# out.write('Cbio_project\tT_CL_BS_ID\tNorm_BS_ID\tFile_Type\tCbio_Tumor_Name\tCbio_Matched_Normal_Name\tFile_Name\n')
# for dpath in ds_list:
#     process_ds(dpath, cav_dict, out)
# out.close()


def process_ds_new(config_data, resourceSet):
    flist = subprocess.check_output('find /Users/kalletlak/Documents/temorary/datasets -name data_clinical_sample.txt', shell=True)

    ds_list = flist.decode().split('\n')
    if ds_list[-1] == '':
        ds_list.pop()

    sbg_api_client = sbg.Api(url='https://cavatica-api.sbgenomics.com/v2', token='296f647e655b4ff2acbc06f92a56b733')
    rsem_files_by_project = {}

    for dpath in ds_list:
        try:
            # last item in find will likely be empty after split
            x = 0
            parts = dpath.split('/')
            cbio_proj = parts[-2]
            # project/disease name should be name of directory hosting datasheet
            sys.stderr.write('Processing ' + cbio_proj + ' project' + '\n')
            cur_ds = open(dpath)
            for i in range(0, 4, 1):
                next(cur_ds)
            ds_head = next(cur_ds)
            ds_header = ds_head.rstrip('\n').split('\t')
            sa_idx = ds_header.index('SAMPLE_ID')
            sp_idx = ds_header.index('SPECIMEN_ID')
            ns_idx = ds_header.index('MATCHED_NORMAL_SPECIMEN_ID')
            studyResources = []
            for entry in cur_ds:
                meta = entry.rstrip('\n').split('\t')
                check = meta[sp_idx].split(';')
                for bs_id in check:
                    # can't tell if RNA or DNA from datasheet, but DNA will be tum + norm bs id, RNA tum + NA
                    key = bs_id + meta[ns_idx]
                    if key in resourceSet:
                        resource = resourceSet[key];
                        resource['SAMPLE_ID'] = meta[sa_idx]
                        studyResources.append(resource)
                    elif bs_id in resourceSet:
                        resource = resourceSet[bs_id];
                        resource['SAMPLE_ID'] = meta[sa_idx]
                        studyResources.append(resource)
                    else:
                        sys.stderr.write('WARNING: ' + key + ' nor ' + bs_id + ' found in data sheet entry\n' + entry + 'Check inputs!\n')

            fileSet = {}
            for studyResource in studyResources:
                for resource in studyResource['resources']:
                    resource.metadata['sample_id'] = studyResource['SAMPLE_ID']
                    name = resource.name
                    name_parts = name.split('.')
                    ext = ".".join(name_parts[1:])

                    if ext in fileSet:
                        fileSet[ext].append(resource)
                    else:
                        fileSet[ext] = [resource]

            maf_extension = config_data['dna_ext_list']['mutation']
            cnv_extension = config_data['dna_ext_list']['copy_number']
            rsem_extension = config_data['rna_ext_list']['expression']
            fusion_extension = config_data['rna_ext_list']['fusion']
            if maf_extension in fileSet:
                add_maf_file(cbio_proj, "/".join(parts[:-1]),  fileSet[maf_extension])
            if cnv_extension in fileSet:
                add_cnv_file(config_data, sbg_api_client, cbio_proj, "/".join(parts[:-1]),  fileSet[cnv_extension])
            if rsem_extension in fileSet:
                rsem_files_by_project[cbio_proj] = fileSet[rsem_extension]
            if fusion_extension in fileSet:
                add_fusion_file(config_data, "/".join(parts[:-1]),  fileSet[fusion_extension])
                # add_fusion_file(cbio_proj, "/".join(parts[:-1]),  fileSet[fusion_extension])

        except Exception as e:
            print(e)
            sys.exit()

    add_rsem_file(config_data, sbg_api_client, rsem_files_by_project)

def add_maf_file(study_name, study_dir, study_maf_files):
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
            process_maf_file('/Users/kalletlak/Documents/temorary/data/' + resource.name, new_maf, t_bs_id, n_bs_id, tid_idx, nid_idx, v_idx, eid_idx, n_idx)
        new_maf.close()
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


def add_cnv_file(config_data, sbg_api_client, study_name, study_dir, study_cnv_resources):
    sys.stderr.write('Processing cnv for ' + study_name + ' project' + '\n')
    data_linear_CNA = open(study_dir + '/data_linear_CNA.txt', 'w')
    data_CNA = open(study_dir + '/data_CNA.txt', 'w')
    cur_cnv_dict = {}
    s_list = []
    high_gain = config_data['cnv_high_gain']

    # TODO: currently if there is a linear cnv data we are adding discrete cna data by default. 
    # We need to check whether *.controlfreec.info.txt files exist in the project

    for resource in study_cnv_resources:
        s_list.append(resource.metadata['sample_id'])
        cur_cnv_dict = process_cnv(config_data, sbg_api_client, study_name, study_dir, resource, cur_cnv_dict)
    data_linear_CNA.write('Hugo_Symbol\t' + '\t'.join(s_list) + '\n')
    data_CNA.write('Hugo_Symbol\t' + '\t'.join(s_list) + '\n')
    for gene in cur_cnv_dict:
        data_linear_CNA.write(gene)
        data_CNA.write(gene)
        for samp in s_list:
            if samp in cur_cnv_dict[gene]:
                data_linear_CNA.write('\t' + str(cur_cnv_dict[gene][samp]['linear_value']))
                if 'discrete_value' in cur_cnv_dict[gene][samp]:
                    data_CNA.write('\t' + str(cur_cnv_dict[gene][samp]['discrete_value']))
                else:
                    data_CNA.write('\t0')
            else:
                data_linear_CNA.write('\t' + str(high_gain))
                data_CNA.write('\t0')
        data_linear_CNA.write('\n')
        data_CNA.write('\n')
    data_linear_CNA.close()
    data_CNA.close()


def process_cnv(config_data, sbg_api_client, study_name, study_dir, resource, cur_cnv_dict):
    try:
        if study_dir[-1] != '/':
            study_dir += '/'
        bedtools = config_data['bedtools']
        cnv_cp_only = config_data['cp_only_script']
        bed_file = config_data['bed_genes']
        hugo_tsv = config_data['hugo_tsv']
        # sys.stderr.write('Processing cna ' + study_name + '\n')
        limit = config_data['cnv_min_len']
        file_name = resource.name
        root_name = os.path.basename(file_name).split('.')[0]
        filtered_bed_file = study_dir + root_name + '.CNVs_' + str(limit) + "_filtered.bed"
        len_fh = open(filtered_bed_file, 'w')
        in_cnv = open('/Users/kalletlak/Documents/temorary/data/'+resource.name)
        skip = next(in_cnv)
        for cnv in in_cnv:
            cnv_data = cnv.rstrip('\n').split('\t')
            if int(cnv_data[2]) - int(cnv_data[1]) >= limit:
                len_fh.write("\t".join(cnv_data[0:5]) + "\n")
        len_fh.close()
        out_file_1 = study_dir + root_name + '.CNVs.Genes'
        out_file = out_file_1 + '.copy_number'
        to_genes_cmd = bedtools + ' intersect -a ' + filtered_bed_file + ' -b ' + bed_file + ' -wb > ' + out_file_1
        subprocess.call(to_genes_cmd, shell=True)
        out_gene_cnv_only = 'perl ' + cnv_cp_only + ' ' + hugo_tsv + ' ' + out_file_1 + ' > ' + out_file
        subprocess.call(out_gene_cnv_only, shell=True)

        info_file_name = file_name.replace("controlfreec.CNVs.p.value.txt", "controlfreec.info.txt")
        response = sbg_api_client.files.query(parent = resource.parent, names=[info_file_name])
        ploidy = None
        if len(response):
            info = response[0].content().split('\n')
            for datum in info:
                if 'Output_Ploidy' in datum:
                    (key, value) = datum.split('\t')
                    ploidy = int(value)
                    break

        high_gain = config_data['cnv_high_gain']
        for entry in open(out_file):
            (gene, gid, linear_value) = entry.rstrip('\n').split('\t')
            if gene not in cur_cnv_dict:
                cur_cnv_dict[gene] = {}
            sample_id = resource.metadata['sample_id']
            if sample_id in cur_cnv_dict[gene]:
                sys.stderr.write('ERROR: Sample ID ' + sample_id + ' already processed.  Back to the drawing board!')
                exit(1)
            else:
                cur_cnv_dict[gene][sample_id] = {}
                if ploidy is not None:
                    linear_value = int(linear_value)
                    discrete_value = linear_value - ploidy
                    min_cn = ploidy * -1
                    # adjust discrete low range only in ploidy > 2
                    if min_cn != -2:
                        discrete_value = np.where(min_cn + 1 <= discrete_value <= -2, -1, discrete_value)
                        discrete_value = np.where(discrete_value == min_cn, -2, discrete_value)
                    discrete_value = np.where(2 <= discrete_value <= high_gain, 1, discrete_value)
                    discrete_value = np.where(discrete_value > high_gain, 2, discrete_value)
                    cur_cnv_dict[gene][sample_id]['discrete_value'] = discrete_value

                cur_cnv_dict[gene][sample_id]['linear_value'] = linear_value
        rm_tmp = 'rm ' + filtered_bed_file + ' ' + out_file_1 + ' ' + out_file
        subprocess.call(rm_tmp, shell=True)
        return cur_cnv_dict
    except Exception as e:
        sys.stderr.write(str(e))
        sys.exit(1)


def mt_collate_df(resource):
    sample_id = resource.metadata['sample_id']
    rsem_file = '/Users/kalletlak/Documents/temorary/data/' + resource.name    
    current = pd.read_csv(rsem_file, sep="\t", index_col=0)
    cur_subset = current.loc[:, ['FPKM']]
    cur_subset.rename(columns={"FPKM": sample_id}, inplace=True)
    return cur_subset


def add_rsem_file(config_data, sbg_api_client, rsem_resources_by_study):
    print('Processing expression files')
    seen_list = []
    df_list = []
    for project in rsem_resources_by_study:
        with concurrent.futures.ThreadPoolExecutor(config_data['cpus']) as executor:
            results = {executor.submit(mt_collate_df, rsem_file): rsem_file for rsem_file in rsem_resources_by_study[project]}
            for result in concurrent.futures.as_completed(results):
                df_list.append(result.result())
    master_tbl = pd.concat(df_list, axis=1)
    del df_list
    # remove duplicate columns/sample_ids
    master_tbl = (master_tbl.loc[:,~master_tbl.columns.duplicated()])
    master_tbl.reset_index(inplace=True)

    gene_id_sym_list = master_tbl['gene_id'].to_list()
    gene_sym_list = []
    for entry in gene_id_sym_list:
        parts = entry.split('_')
        gene_sym = '_'.join(parts[1:])
        gene_sym_list.append(gene_sym)
    master_tbl['Hugo_Symbol'] = gene_sym_list
    master_tbl.drop(columns =["gene_id"], inplace = True)
    dup_hugo_tbl = master_tbl[master_tbl.duplicated(['Hugo_Symbol'])]
    rpt_sym = dup_hugo_tbl.Hugo_Symbol.unique()

    sample_list = (master_tbl.columns.difference(['Hugo_Symbol']))
    for sym in rpt_sym:
        gene_eval = master_tbl.loc[master_tbl['Hugo_Symbol'] == sym]
        mean_expr = gene_eval[sample_list].mean(axis=1).sort_values(ascending=False)
        to_drop = list(mean_expr.index)[1:]
        master_tbl.drop(to_drop, inplace=True)
    
    master_tbl.set_index('Hugo_Symbol', inplace=True)
    gene_sym_list = master_tbl.index

    sys.stderr.flush()

    for project in rsem_resources_by_study:
        sampleIds = []
        for resource in rsem_resources_by_study[project]:
            sampleIds.append(resource.metadata['sample_id'])
        expr_fname = '/Users/kalletlak/Documents/temorary/datasets/' + project + '/data_rna_seq_v2_mrna.txt'
        master_tbl[sampleIds].to_csv(expr_fname, sep='\t', mode='w', index=True)


    z_scored = stats.zscore(np.log2(np.array(master_tbl + 1)), axis = 1)
    del master_tbl
    master_zscore_log = pd.DataFrame(z_scored, index=gene_sym_list, columns=sample_list)
    # this may be memory-intensive for some insane reson...
    sys.stderr.write("Replacing NaN with 0\n")
    sys.stderr.flush()
    master_zscore_log.fillna(0, inplace=True)
    sys.stderr.write('Outputing z scored results\n')
    sys.stderr.flush()

    for project in rsem_resources_by_study:
        sampleIds = []
        for resource in rsem_resources_by_study[project]:
            sampleIds.append(resource.metadata['sample_id'])
        zscore_fname = '/Users/kalletlak/Documents/temorary/datasets/' + project + '/data_rna_seq_v2_mrna_median_Zscores.txt'
        master_zscore_log[sampleIds].to_csv(zscore_fname, sep='\t', mode='w', index=True)


# def add_fusion_file(study_name, study_dir, study_cnv_resources):
#     sys.stderr.write('Processing fusion for ' + study_name + ' project' + '\n')
#     if study_dir[-1] != '/':
#             study_dir += '/'
#     frame_list = []
#     for resource in study_cnv_resources:
#         file_name = resource.name
#         ann_file = pd.read_csv('/Users/kalletlak/Documents/temorary/data/'+resource.name, sep="\t")
#         ann_file['Tumor_Sample_Barcode'] = resource.metadata['sample_id']
#         frame_list.append(ann_file)
#     concat_frame = pd.concat(frame_list)
#     del frame_list
#     fusion_data = concat_frame[['Gene1A', 'Tumor_Sample_Barcode', 'FusionName', 'Caller', 'Fusion_Type']]
#     del concat_frame
#     fusion_data = fusion_data.groupby(['Gene1A', 'Tumor_Sample_Barcode', 'FusionName', 'Fusion_Type'])['Caller'].unique().apply(', '.join).reset_index()
#     fusion_data['Caller'] = fusion_data['Caller'].str.upper()
#     fusion_data['Center'] = 'Sequencing center'
#     fusion_data.rename(columns={"Gene1A": "Hugo_Symbol", "FusionName": "Fusion", "Caller": "Method",
#     "Fusion_Type": "Frame"}, inplace=True)
#     fusion_data['DNA_support'] = "no"
#     fusion_data['RNA_support'] = "yes"
#     # create blank entrez ID column so that 3' gene names can be searched
#     fusion_data['Entrez_Gene_Id'] = ""
#     # Reorder table
#     order_list = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'Tumor_Sample_Barcode', 'Fusion','DNA_support', 'RNA_support', 'Method', 'Frame']
#     fusion_data = fusion_data[order_list]
#     fusion_data.set_index('Hugo_Symbol', inplace=True)

#     fus_file = open(study_dir + 'data_fusions.txt', 'w')
#     fusion_data.reset_index(inplace=True)
#     fus_head = "\t".join(fusion_data.columns) + "\n"
#     fus_file.write(fus_head)
#     fus_data = list(fusion_data.values)
#     # "hack" to allow 3' end of fusion to be searched
#     for data in fus_data:
#         fus_file.write("\t".join(data) + "\n")
#         (a, b) = data[4].split('--')
#         if a != b:
#             data[0] = b
#             fus_file.write("\t".join(data) + "\n")
#     fus_file.close()