import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Script to pick out samples from CNVS file when sample type in sample map file is DNA")

    parser.add_argument('-ID','--sample_map',help="File to read in sample ID")
    parser.add_argument('-cnv','--cnvs',help="File to pick cnv samples")
    parser.add_argument('-o','--output_file',help="Name of output file")

    args=parser.parse_args()

    sample_map = pd.read_csv(args.ID, sep='\t')

    #predicted = pd.read_csv("openpbta.predicted_cnv.txt", sep='\t')
    predicted = pd.read_csv(args.cnvs, sep='\t')
    cbio_ID_DNA=sample_map.loc[sample_map['Sample Type']=='DNA']['Cbio ID']

    #print(len(cbio_ID_DNA))
    #print(predicted.shape)
    add=[]
    a=0
    for col in cbio_ID_DNA:
        a=0
        for row in predicted.columns:
            if row == col:
                a=1
        if a!=1:
            add.append(col)
            predicted[col]=0

    #print(predicted.shape)
    with open(args.output_file,'w') as write_tsv:
        write_tsv.write(predicted.to_csv(sep='\t', index=False))    
