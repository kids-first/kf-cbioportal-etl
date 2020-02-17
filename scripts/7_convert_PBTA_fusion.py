import sys
import argparse
import json
import os
import pandas as pd
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert openPBTA fusion table to cbio format.')
    parser.add_argument('-t', '--table', action='store', dest='table',
                    help='Table with cbio project, kf bs ids, cbio IDs, and file names')
    parser.add_argument('-f', '--fusion-file', action='store', dest='fusion_file', help='openPBTA fusion file')
    parser.add_argument('-s', '--center-file', action='store', dest='sq_file', help='File with BS IDs and sequencing centers')
    parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types '
                                                                                   'and data locations')
    args = parser.parse_args()

    with open(args.config_file) as f:
        config_data = json.load(f)

    out_dir = 'merged_fusion/'
    try:
        os.mkdir(out_dir)
    except:
        sys.stderr.write('output dir already exists\n')
        sys.stderr.flush()

    fusion_data = pd.read_csv(args.fusion_file, sep="\t")
    sq_info = pd.read_csv(args.sq_file, sep="\t")
    all_file_meta = pd.read_csv(args.table, sep="\t")
    # entrez_data = pd.read_csv(config_data['entrez_tsv'], sep="\t")

    # deal only with RNA metadata
    rna_subset = all_file_meta.loc[all_file_meta['File_Type'] == 'rsem']
    rna_subset.rename(columns={"T_CL_BS_ID": "Sample"}, inplace=True)
    rna_subset.set_index('Sample', inplace=True)

    # drop phosho entrez IDs
    # no_phospho = entrez_data[entrez_data['ENTREZ_GENE_ID'] > 0]
    # Drop duplicates
    # no_phospho.drop_duplicates(subset ="HUGO_GENE_SYMBOL", keep = False, inplace = True)
    # no_phospho.rename(columns={"HUGO_GENE_SYMBOL": "Hugo_Symbol"}, inplace=True)
    # no_phospho.set_index('Hugo_Symbol', inplace=True)
    # del entrez_data
    # Get releavent columns
    cbio_master = fusion_data[['Gene1A', 'Sample', 'FusionName', 'CalledBy', 'Fusion_Type']]
    cbio_master.set_index('Sample', inplace=True)
    sq_info.rename(columns={"BS_ID": "Sample"}, inplace=True)
    sq_info.set_index('Sample', inplace=True)
    # get cbio IDs and seq center
    cbio_master = pd.merge(cbio_master, rna_subset[['Cbio_Tumor_Name']], on='Sample', how='inner')
    cbio_master = pd.merge(cbio_master, sq_info[['SQ_Value']], on='Sample', how='inner')
    # drop bs ids
    cbio_master.reset_index(inplace=True)
    cbio_master.drop('Sample', axis=1, inplace=True)
    cbio_master.rename(columns={"Gene1A": "Hugo_Symbol", "FusionName": "Fusion", "CalledBy": "Method",
    "Fusion_Type": "Frame", "Cbio_Tumor_Name": "Tumor_Sample_Barcode", "SQ_Value": "Center"}, inplace=True)
    cbio_master['DNA_support'] = "no"
    cbio_master['RNA_support'] = "yes"
    # cbio_master = pd.merge(cbio_master, no_phospho, on='Hugo_Symbol', how='left')
    cbio_master = cbio_master.replace(np.nan, '', regex=True)
    cbio_master.rename(columns={"ENTREZ_GENE_ID": "Entrez_Gene_Id"}, inplace=True)
    # Reorder table
    cbio_master['Entrez_Gene_Id'] = ""
    order_list = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'Tumor_Sample_Barcode', 'Fusion','DNA_support', 'RNA_support', 'Method', 'Frame']
    cbio_master = cbio_master[order_list]
    project_list = rna_subset.Cbio_project.unique()
    cbio_master.set_index('Hugo_Symbol', inplace=True)
    # cbio_master[['Fusion']] = cbio_master[['Fusion']].replace(to_replace=r'--', value='-', regex=True)
    # remove trailing .0 from entrez ID being treated as number
    # cbio_master['Entrez_Gene_Id']=cbio_master.Entrez_Gene_Id.apply(lambda x: str(x))
    # cbio_master['Entrez_Gene_Id'] = cbio_master['Entrez_Gene_Id'].str.replace('.0$', '')
    for project in project_list:
        sub_sample_list = list(rna_subset.loc[rna_subset['Cbio_project'] == project, 'Cbio_Tumor_Name'])
        fus_fname = out_dir + project + '.fusions.txt'
        fus_file = open(fus_fname, 'w')
        fus_tbl = cbio_master[cbio_master.Tumor_Sample_Barcode.isin(sub_sample_list)]
        fus_tbl.reset_index(inplace=True)
        fus_head = "\t".join(fus_tbl.columns) + "\n"
        fus_file.write(fus_head)
        fus_data = list(fus_tbl.values)
        # "hack" to allow 3' end of fusion to be searched
        for data in fus_data:
            fus_file.write("\t".join(data) + "\n")
            (a, b) = data[4].split('--')
            if a != b:
                data[0] = b
                fus_file.write("\t".join(data) + "\n")
        fus_file.close()
        # cbio_master[cbio_master.Tumor_Sample_Barcode.isin(sub_sample_list)].to_csv(fus_fname, sep='\t', mode='w', index=True)
