#!/usr/bin/env python3
import sys
import pandas as pd
from .meta_file_utils import add_meta_file


def add_maf_file(config_data, study_name, study_dir, study_maf_files):
    maf_exc =["Silent", "Intron", "IGR", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "RNA"]
    maf_df_list = []
    try:
        sys.stderr.write('Processing for maf ' + study_name + ' project' + '\n')
        for resource in study_maf_files:
            maf_df = pd.read_csv(config_data['data_files'] + resource.name, sep="\t", skiprows=1)
            maf_df['Tumor_Sample_Barcode'] = resource.metadata['Kids First Biospecimen ID Tumor']
            maf_df['Matched_Norm_Sample_Barcode'] = resource.metadata['Kids First Biospecimen ID Normal']
            maf_df = maf_df[~maf_df['Variant_Classification'].isin(maf_exc)]
            maf_df_list.append(maf_df)

        maf_data = pd.concat(maf_df_list)
        del maf_data['Entrez_Gene_Id'] # is this required?
        maf_data.to_csv(study_dir + "/data_mutations_extended.txt", sep="\t", index=False)

        #Add meta file
        add_meta_file(study_name, config_data["metadata"]["mutation"]["meta_attr"], study_dir + '/meta_mutations_extended.txt')

    except Exception as e:
        print(e)
        sys.exit()