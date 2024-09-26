#!/usr/bin/env python3
"""
Likely a one-off script. Compare CN matrices from a previous load to a new one, summarize difference, and output number of hotspot changes
Usage:
  python3 compare_cn_files.py old.tsv new.tsv hotspots.txt
"""

import pandas as pd
import sys
import numpy as np

old = pd.read_csv(sys.argv[1], sep="\t")
old.set_index("Hugo_Symbol", inplace=True)
old.sort_index(inplace=True)

new = pd.read_csv(sys.argv[2], sep="\t")
new.set_index("Hugo_Symbol", inplace=True)
new.sort_index(inplace=True)

ne_stacked = (old != new).stack()
changed = ne_stacked[ne_stacked]

difference_locations = np.where(old != new)

changed_from = old.values[difference_locations]

changed_to = new.values[difference_locations]
result =  pd.DataFrame({'from': changed_from, 'to': changed_to}, index=changed.index)

with open(sys.argv[3]) as f:
    hotspots = f.read().splitlines()

print("HotSpot\tCount Added\tCount Dropped")
new_added = result[result['from'] == 0]
print("{} entries will be added".format(len(new_added)), file=sys.stderr)
new_rm = result[result['to'] == 0]
print("{} entries will be dropped".format(len(new_rm)), file=sys.stderr)
clean_list = []
for gene in hotspots:
    add = 0
    drop = 0
    try:
        add = len(new_added.loc[[gene]])
    except Exception as e:
        print(e, file=sys.stderr)
    try:
        drop = len(new_rm.loc[[gene]])
    except Exception as e:
        print(e, file=sys.stderr)
    if add or drop:
        print("{}\t{}\t{}".format(gene, add, drop))
        clean_list.append(gene)
