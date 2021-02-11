#!/usr/bin/env python3


def add_meta_file(cancer_study_identifier, meta_information, meta_file_path):
    meta_file = open(meta_file_path, 'w')
    meta_file.write("cancer_study_identifier: " + cancer_study_identifier + "\n")
    for mkey in meta_information:
        meta_file.write(mkey + ": " + meta_information[mkey] + "\n")
    meta_file.close()
