import sys
import argparse
import os
import pandas as pd
import re
import pdb

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate extra cases lists based on cBio datasheet')
    parser.add_argument('-d', '--datasheet', action='store', dest='datasheet',
                    help='cBio sample datasheet')
    parser.add_argument('-s', '--study-id', action='store', dest='study_id',
                    help='cBio cancer_study_identifier')
    
    args = parser.parse_args()
    # skip first 4 header rows
    cbio_ds = pd.read_csv(args.datasheet, sep="\t", skiprows=4)
    # iterate through histologies with at least 5 samples
    hist_list = cbio_ds.HISTOLOGY.unique()
    cohort_list = cbio_ds.COHORT.unique()
    for hist in hist_list:
        samp_ids = cbio_ds[cbio_ds['HISTOLOGY'].isin([hist])].SAMPLE_ID.to_list()
        if len(samp_ids) >= 5:
            # print(hist)
            try:
                suffix = hist.lower().replace(" ", "_").replace("/", "_")
                cur = open("cases_" + suffix + ".txt", "w")
                cur.write("cancer_study_identifier: " + args.study_id + "\n" 
                + "stable_id: " + args.study_id + "_" + suffix + "\n" 
                + "case_list_name: " + hist + " samples\n" 
                + "case_list_description:  " + hist + " samples (" + str(len(samp_ids)) + ")\n"
                + "case_list_category: other\ncase_list_ids: " + "\t".join(samp_ids) + "\n")
                cur.close()
                # Also collapse on the cohort level
                for cohort in cohort_list:
                    subset = cbio_ds.loc[cbio_ds['COHORT'] == cohort,]
                    subset_samp = subset[subset['HISTOLOGY'].isin([hist])].SAMPLE_ID.to_list()
                    if len(subset_samp) >= 5:
                        cohort_stable = cohort.lower().replace(" ", "_").replace("/", "_")
                        cur = open("cases_" + cohort_stable + "_" + suffix + ".txt", "w")
                        cur.write("cancer_study_identifier: " + args.study_id + "\n" 
                        + "stable_id: " + args.study_id + "_" + cohort_stable +  "_" + suffix + "\n" 
                        + "case_list_name: " + cohort + " " + hist + " samples\n" 
                        + "case_list_description:  " + cohort + " " + hist + " samples (" + str(len(subset_samp)) + ")\n"
                        + "case_list_category: other\ncase_list_ids: " + "\t".join(subset_samp) + "\n")
                        cur.close()
                    else:
                        sys.stderr.write("Skipping " + hist + " for cohort " + cohort + ", less than 5 entries\n")

            except Exception as e:
                print(str(e))
                sys.stderr.write("Error encountered - likely a blank histology, skipping\n")
        else:
            sys.stderr.write("Skipping " + hist + " less than 5 entries\n")
    #pdb.set_trace()
    #hold=1