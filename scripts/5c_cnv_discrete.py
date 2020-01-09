#!/usr/bin/env python3

import sevenbridges as sbg
import sys
import argparse
from sevenbridges.http.error_handlers import rate_limit_sleeper, maintenance_sleeper
import concurrent.futures
import json
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Convert merged cnv values to discrete coded values.')
parser.add_argument('-f', '--merged_cnv_file', action='store', dest='cnv_file', help='merged cnv file')
parser.add_argument('-j', '--config', action='store', dest='config_file', help='json config file with data types and '
                                                                               'data locations')
parser.add_argument('-m', '--info_manifest', action='store', dest='manifest', help='cavatica cfree info file manifest')
parser.add_argument('-p', '--profile', action='store', dest='profile', help='cavatica profile name. requires '
                                                                            '.sevenbridges/credentials file be present')


def mt_adjust_cn(obj):
    try:
        samp_id = obj.resource.metadata['case_id']
        info = obj.resource.content().split('\n')
        info_dict = {}
        for datum in info:
            try:
                (key, value) = datum.split('\t')
                info_dict[key] = value
            except Exception as e:
                sys.stderr.write(str(e) + " entry could not be split. skipping!\n")
        ploidy = int(info_dict['Output_Ploidy'])
        data[samp_id] = data[samp_id] - ploidy
        min_cn = ploidy * -1
        # adjust discrete low range only in ploidy > 2
        if min_cn != -2:
            data[samp_id] = np.where(data[samp_id].between((min_cn + 1),-2), -1, data[samp_id])
            data[samp_id] = np.where(data[samp_id] == min_cn, -2, data[samp_id])
        data[samp_id] = np.where(data[samp_id].between(2, high_gain), 1, data[samp_id])
        data[samp_id] = np.where(data[samp_id] > high_gain, 2, data[samp_id])
        return [0, obj.resource.id]

    except Exception as E:
        sys.stderr.write(str(E) + "\n")
        sys.stderr.write('Error processing sample entry for ' + obj.resource.id + '\n')
        return [1, obj.resource.id]


args = parser.parse_args()
with open(args.config_file) as f:
    config_data = json.load(f)

config = sbg.Config(profile=args.profile)
api = sbg.Api(config=config, error_handlers=[
                sbg.http.error_handlers.rate_limit_sleeper,
                sbg.http.error_handlers.maintenance_sleeper])

data = pd.read_csv(args.cnv_file, sep="\t")
data.set_index('Hugo_Symbol')

samp_list = list(data.columns)[1:]
manifest = pd.read_csv(args.manifest)
# may need to change this for other projects :(
manifest.set_index('case_id', inplace=True)
bulk_ids = []
#fid_dict = {}
for samp_id in samp_list:
    bulk_ids.append(manifest.loc[samp_id]['id'])
    # fid_dict[samp_id] = manifest.loc[samp_id]['id']
file_objs = []
max_j = 100
total = len(bulk_ids)
for i in range(0, total, max_j):
    uset = i + max_j
    sys.stderr.write('Processing ' + str(uset) + ' set\n')
    if uset > total:
        uset = total
    subset = api.files.bulk_get(files=bulk_ids[i:uset])
    file_objs.extend(subset)
# may want to make a config file option
high_gain = 4

x = 1
m = 100
with concurrent.futures.ThreadPoolExecutor(config_data['threads']) as executor:
    results = {executor.submit(mt_adjust_cn, obj): obj for obj in file_objs}
    for result in concurrent.futures.as_completed(results):
        if result.result()[0] == 1:
            'Had trouble processing object ' + result.result([1] + '\n')
            sys.exit(1)
        if x % m == 0:
            sys.stderr.write('Processed ' + str(x) + ' samples\n')
            sys.stderr.flush()
        x += 1
sys.stderr.write('Conversion completed.  Writing results to file\n')
data.to_csv(r'discrete_cnvs.txt', sep='\t', mode='w', index=False)
