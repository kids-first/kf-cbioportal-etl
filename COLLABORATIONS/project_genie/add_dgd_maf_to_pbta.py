import sys
import os
import argparse
import pdb

parser = argparse.ArgumentParser(description='Output fields from maf file based on header - meant to be appended to an existing file!')
parser.add_argument('-i', '--header', action='store', dest='header', help='File with maf header only')
parser.add_argument('-m', '--maf-dir', action='store', dest='maf_dir', help='maf file directory')

args = parser.parse_args()
maf_list = os.listdir(args.maf_dir)
header_file = open(args.header)
skip = next(header_file)
head = next(header_file)
# store header fields as list, will use to get position in new file, if exists
# will output blank if not in new file
header = head.rstrip('\n').split('\t')
h_dict = {}
for item in header:
    h_dict[item] = None
for maf in maf_list:
    sys.stderr.write("Processing " + maf + "\n")
    cur = open(args.maf_dir + "/" + maf)
    skip = next(cur)
    m_head = next(cur)
    m_header = m_head.rstrip('\n').split('\t')
    for i in range(len(m_header)):
        if m_header[i] in h_dict:
            h_dict[m_header[i]] = i
    # only print items in original headr and in same order, else print blank
    for data in cur:
        to_print = []
        datum = data.rstrip('\n').split('\t')
        for item in header:
            if h_dict[item] != None:
                to_print.append(datum[h_dict[item]])
            else:
                to_print.append("")
        sys.stdout.write("\t".join(to_print) + "\n")
    sys.stderr.write("Processed " + maf + "\n")