#!/usr/bin/env python3
import sys
import pandas as pd
from requests import request

from .meta_file_utils import add_meta_file

def get_sequencing_center_info(config_data, bs_id):
    bs_url =  config_data['kf_url'] + '/biospecimens/' + bs_id
    bs_info = request('GET', bs_url)
    bs_info.json()
    if bs_info.json()['_status']['code'] == 404:
        sys.stderr.write(bs_id + ' not found!\n')
        return ''
    sequencing_center_url = config_data['kf_url'] + bs_info.json()['_links']['sequencing_center']
    sequencing_center_info = request('GET', sequencing_center_url)
    return sequencing_center_info.json()['results']['name']


def add_fusion_file(config_data, study_dir, study_resources):
    if study_dir[-1] != '/':
        study_dir += '/'
    
    study_name = study_dir.split('/')[-2]
    
    sys.stderr.write('Processing fusion for ' + study_name + ' project' + '\n')
    data_files_directory = config_data['data_files']
    frame_list = []
    for resource in study_resources:
        ann_file = pd.read_csv(data_files_directory+resource.name, sep="\t")
        ann_file['Tumor_Sample_Barcode'] = resource.metadata['sample_id']
        ann_file['Center'] = get_sequencing_center_info(config_data, ann_file.at[0, 'Sample'])
        frame_list.append(ann_file)
    fusion_data = pd.concat(frame_list)
    del frame_list
    fusion_data = fusion_data[['Gene1A', 'Tumor_Sample_Barcode', 'FusionName', 'Caller', 'Fusion_Type', 'Center']]
    fusion_data = fusion_data.groupby(['Gene1A', 'Tumor_Sample_Barcode', 'FusionName', 'Fusion_Type', 'Center'])['Caller'].unique().apply(', '.join).reset_index()
    fusion_data['Caller'] = fusion_data['Caller'].str.upper()
    fusion_data.rename(columns={"Gene1A": "Hugo_Symbol", "FusionName": "Fusion", "Caller": "Method", "Fusion_Type": "Frame"}, inplace=True)
    fusion_data['DNA_support'] = "no"
    fusion_data['RNA_support'] = "yes"
    # create blank entrez ID column so that 3' gene names can be searched
    fusion_data['Entrez_Gene_Id'] = ""
    # Reorder table
    order_list = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'Tumor_Sample_Barcode', 'Fusion','DNA_support', 'RNA_support', 'Method', 'Frame']
    fusion_data = fusion_data[order_list]
    fusion_data.set_index('Hugo_Symbol', inplace=True)

    #Add data file
    fus_file = open(study_dir + 'data_fusions.txt', 'w')
    fusion_data.reset_index(inplace=True)
    fus_head = "\t".join(fusion_data.columns) + "\n"
    fus_file.write(fus_head)
    #fus_data = list(fusion_data.values)
    # "hack" to allow 3' end of fusion to be searched
    for data in fusion_data.values.tolist():
        fus_file.write("\t".join(data) + "\n")
        (a, b) = data[4].split('--')
        if a != b:
            data[0] = b
            fus_file.write("\t".join(data) + "\n")
    fus_file.close()

    #Add meta file
    add_meta_file(study_name, config_data["metadata"]["fusion"]["meta_attr"], study_dir + 'meta_fusions.txt')