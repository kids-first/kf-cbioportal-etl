import sys
import argparse
import pandas as pd
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add fusion flag to data clinical sample')
    parser.add_argument('-f', '--fusion-dir', action='store', dest='fusion_dir',
                    help='Dir with fusions')
    parser.add_argument('-d', '--datasheet-list', action='store', dest='datasheet_list', help='cbio data sample sheet list')
    parser.add_argument('-o', '--out-dir', action='store', dest='out_dir', help='output dir - don\'t make it the same as current datasheet location!!!')

    args = parser.parse_args()
    added_header = ['WITH_FUSION_DATA',
    'Indicator that sample has fusions associated with it',
    'STRING',
    '1',
    'WITH_FUSION_DATA']
    os.mkdir(args.out_dir)
    datasheet_list = open(args.datasheet_list)
    for fname in datasheet_list:
        datasheet = open(fname.rstrip('\n'))
        parts = fname.rstrip('\n').split("/")
        fusion_data = pd.read_csv(args.fusion_dir + "/" + parts[-2] + ".fusions.txt", sep="\t")
        sample_list = fusion_data.Tumor_Sample_Barcode.unique()
        cur_dir = args.out_dir + "/" + parts[-2]
        os.mkdir(cur_dir)
        out = open(cur_dir + "/" + parts[-1], "w")
        for i in range(0, len(added_header)):
            head = next(datasheet)
            out.write(head.rstrip('\n') + '\t' + added_header[i] + '\n')
        for line in datasheet:
            out.write(line.rstrip('\n'))
            data = line.split('\t')
            if data[1] in sample_list:
                out.write('\tYES\n')
            else:
                out.write('\tNO\n')
        out.close()