#!/usr/bin/env python3
import sys
import os
import subprocess
import numpy as np

from .file_utils import write_meta_file

def add_cnv_file(config_data, sbg_api_client, study_name, study_dir, study_cnv_resources):
    sys.stderr.write('Processing cnv for ' + study_name + ' project' + '\n')
    cur_cnv_dict = {}
    s_list = []
    high_gain = config_data['cnv_high_gain']

    for resource in study_cnv_resources:
        s_list.append(resource.metadata['sample_id'])
        cur_cnv_dict = process_cnv(config_data, sbg_api_client, study_name, study_dir, resource, cur_cnv_dict)

    data_linear_CNA = open(study_dir + '/data_linear_CNA.txt', 'w')
    data_CNA = open(study_dir + '/data_CNA.txt', 'w')
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

    #Add meta file
    write_meta_file(study_name, config_data["metadata"]["cnv"]["linear"]["meta_attr"], study_dir + '/meta_linear_CNA.txt')
    write_meta_file(study_name, config_data["metadata"]["cnv"]["discrete"]["meta_attr"], study_dir + '/meta_CNA.txt')


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
        in_cnv = open(config_data['data_files']+resource.name)
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

        # aasuming *.controlfreec.info.txt exists for every *.controlfreec.CNVs.p.value.txt
        info_file_name = file_name.replace("controlfreec.CNVs.p.value.txt", "controlfreec.info.txt")
        response = sbg_api_client.files.query(parent = resource.parent, names=[info_file_name])
        # TODO: revert back to None once data files are updated
        ploidy = config_data['ploidy']
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
                sys.exit(1)
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
        sys.stderr.write('Error ' + str(e) + ' occurred while trying to process cnv file')
        sys.exit(1)