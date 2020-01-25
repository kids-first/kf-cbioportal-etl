import sys
import argparse
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add fusion flag to data clinical sample')
    parser.add_argument('-f', '--fusion-file', action='store', dest='fusion_file',
                    help='File with fusions')
    parser.add_argument('-d', '--datasheet', action='store', dest='datasheet', help='cbio data sample sheet')

    args = parser.parse_args()
    fusion_data = pd.read_csv(args.fusion_file, sep="\t")
    added_header = ['With Fusion Data',
    'Indicator that sample has fusions associated with it',
    'STRING',
    '1',
    'WITH_FUSION_DATA']

    sample_list = fusion_data.Tumor_Sample_Barcode.unique()
    datasheet = open(args.datasheet)
    for i in range(0, len(added_header)):
        head = next(datasheet)
        sys.stdout.write(head.rstrip('\n') + '\t' + added_header[i] + '\n')
    for line in datasheet:
        sys.stdout.write(line.rstrip('\n'))
        data = line.split('\t')
        if data[1] in sample_list:
            sys.stdout.write('\tYES\n')
        else:
            sys.stdout.write('\tNO\n')