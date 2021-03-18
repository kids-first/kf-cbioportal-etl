#!/usr/bin/env python3


def write_meta_file(cancer_study_identifier, meta_information, meta_file_path):
    meta_file = open(meta_file_path, 'w')
    meta_file.write("cancer_study_identifier: " + cancer_study_identifier + "\n")
    for mkey in meta_information:
        meta_file.write(mkey + ": " + meta_information[mkey] + "\n")
    meta_file.close()


def write_case_list(case_file, attr_dict, canc_study_id, sample_list):
    case_file = open(case_file, "w")
    case_file.write("cancer_study_identifier: " + canc_study_id + "\n")
    key_list = list(attr_dict)
    case_file.write(key_list[0] + ": " + canc_study_id + "_" + attr_dict[key_list[0]] + "\n")
    for i in range(1, len(key_list), 1):
        to_write = key_list[i] + ": " + attr_dict[key_list[i]]
        if key_list[i] == "case_list_description":
            to_write += " (" + str(len(sample_list)) + ")"
        case_file.write(to_write + "\n")
    case_file.write("case_list_ids: " + "\t".join(sample_list) + "\n")
    case_file.close()
