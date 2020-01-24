import sys
import argparse
import json
import pandas as pd
import numpy as np
import pdb

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

    pdb.set_trace()
    fusion_data = pd.read_csv(args.fusion_file, sep="\t")
    sq_info = pd.read_csv(args.sq_file, sep="\t")
    all_file_meta = pd.read_csv(args.table, sep="\t")
    entrez_data = pd.read_csv(config_data['hugo_tsv'], sep="\t")

    # deal only with RNA metadata
    rna_subset = all_file_meta.loc[all_file_meta['File_Type'] == 'rsem'][['T_CL_BS_ID', 'Cbio_Tumor_Name']]
    rna_subset.set_index('T_CL_BS_ID', inplace=True)

    # drop phosho entrez IDs
    no_phospho = entrez_data[entrez_data['ENTREZ_GENE_ID'] < 0]
    # Drop duplicates
    no_phospho.drop_duplicates(subset ="HUGO_GENE_SYMBOL", keep = False, inplace = True)
    no_phospho.set_index('HUGO_GENE_SYMBOL', inplace=True)
    del entrez_data
    # Get releavent columns
    cbio_master = fusion_data[['Gene1A', 'Sample', 'FusionName', 'CalledBy', 'Fusion_Type']]
    cbio_master.set_index('Sample', inplace=True)
    sq_info.set_index('BS_ID', inplace=True)
    # get cbio IDs and seq center
    cbio_master = pd.concat([cbio_master, rna_subset], axis=1)
    cbio_master = pd.concat([cbio_master, sq_info[['SQ_Value']]], axis=1)
    # drop bs ids
    cbio_master.reset_index(inplace=True)
    cbio_master.drop('Sample', axis=1, inplace=True)
    cbio_master.rename(columns={"Gene1A": "Hugo_Symbol", "FusionName": "Fusion", "CalledBy": "Method",
    "Fusion_Type": "Frame", "Cbio_Tumor_Name": "Tumor_Sample_Barcode", "SQ_Value", "Center"}, inplace=True)
    cbio_master['DNA_support'] = "no"
    cbio_master['RNA_support'] = "yes"
    gene_list = cbio_master.Hugo_Symbol.uniq()
    





