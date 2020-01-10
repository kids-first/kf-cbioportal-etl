import sys
import argparse
import os
import concurrent.futures
import json
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description='Merge rsem files using cavatica file info.')
parser.add_argument('-f', '--merged_rsem', action='store', dest='merged_rsem',
                    help='File from step 6a, merged rsem file')
parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types '
                                                                                'and data locations')
args = parser.parse_args()

with open(args.config_file) as f:
    config_data = json.load(f)

data = pd.read_csv(args.merged_rsem, sep="\t")
data.set_index('Hugo_Symbol')