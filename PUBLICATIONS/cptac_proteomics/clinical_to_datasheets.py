import sys
import argparse
import json
import math


def format_desc(entry):
    desc = entry.replace(".", " ")
    return desc

def format_colname(entry):
    col_name = entry.replace(".", "_").upper()
    return col_name


def init_list(header_list, entry):
    header_list.append([])
    header_list[0].append(entry)
    header_list.append([])
    header_list[1].append(format_desc(entry))
    header_list.append([])
    header_list[2].append(h_dict[entry][3])
    header_list.append([])
    header_list[3].append("1")
    header_list.append([])
    if h_dict[entry][0] == "":
        header_list[4].append(format_colname(entry))
    else:
        header_list[4].append(h_dict[entry][0])
    return header_list


def build_header(header_list, entry):

    if not header_list:
        header_list = init_list(header_list, entry)
    else:
        header_list[0].append(entry)
        header_list[1].append(format_desc(entry))
        header_list[2].append(h_dict[entry][3])
        header_list[3].append("1")
        if h_dict[entry][0] == "":
            header_list[4].append(format_colname(entry))
        else:
            header_list[4].append(h_dict[entry][0])

    return header_list

parser = argparse.ArgumentParser(description='Script to convert clinical data to cbio clinical data sheets')
parser.add_argument('-f', '--header-file', action='store', dest="head" ,help='tsv file with input file original sample '
                    'names, output sheet flag, and conversion')
parser.add_argument('-c', '--clinical-data', action='store', dest='clin',
                    help='Input clinical data sheet')

args = parser.parse_args()

# read in header key file to assign some of the clinical data headers to cbio standard ones and to all a data type and sheet location
# File header is: original	MAP TO	sample	patient	type
h_dict = {}
h_file = open(args.head)
skip = next(h_file)
for line in h_file:
    info = line.rstrip('\n').split('\t')
    h_dict[info[0]] = info[1:]

clin_data = open(args.clin)
head = next(clin_data)
"""
Data sheet header format:
["#Patient Identifier\tSample Identifier\tSPECIMEN_ID\tCANCER_TYPE\tCANCER_TYPE_DETAILED\tTUMOR_TISSUE_SITE\tTUMOR_TYPE\tSAMPLE_TYPE\tMATCHED_NORMAL_SAMPLE_ID\tMATCHED_NORMAL_SPECIMEN_ID\tCBTTC_PAIRED_IDS",
    "#Patient identifier\tSample Identifier using external_sample_id\tkfdrc tumor biopsecimen ID\tStudy-defined cancer type\tStudy-defined cancer type detail\ttumor tissue location\tprimary v metastatic tumor designation\tpatient tissue sample or cell line\tmatched normal external_sample_id\tkfdrc matched normal biospecimen ID\tOriginal CBTTC DNA/RNA pair ids, where applicable",
    "#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING",
    "#1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1"
    "PATIENT_ID\tSAMPLE_ID\tSPECIMEN_ID\tCANCER_TYPE\tCANCER_TYPE_DETAILED\tTUMOR_TISSUE_SITE\tTUMOR_TYPE\tSAMPLE_TYPE\tMATCHED_NORMAL_SAMPLE_ID\tMATCHED_NORMAL_SPECIMEN_ID\tCBTTC_PAIRED_IDS\n"],

"""
header = head.rstrip('\n').split('\t')
# name of columns to use to calculate AGE value
age_names = ["Overall.Survival", "Age.at.Last.Known.Clinical.Status"]
outs = {}
sample_out = open("data_clinical_sample.txt", "w")
patient_out = open("data_clinical_patient.txt", "w")
# initialize first entry, as it has all leading "#""
patient_head_list = []
sample_head_list = []
for entry in header:
    if h_dict[entry][1] == '1':
        sample_head_list = build_header(sample_head_list, entry)
    if h_dict[entry][2] == '1':
        patient_head_list = build_header(patient_head_list, entry)
# add at end specially calculated age column
age_header = ["AGE", "Age at which the condition or disease was first diagnosed, in years", "NUMBER", "1", "AGE"]
for i in range(len(age_header)):
    patient_head_list[i].append(age_header[i])
# output headers to files
for i in range(0, len(age_header)-1, 1):
    sample_out.write("#" + "\t".join(sample_head_list[i]) + "\n")
    patient_out.write("#" + "\t".join(patient_head_list[i]) + "\n")
sample_out.write("\t".join(sample_head_list[-1]) + "\n")
patient_out.write("\t".join(patient_head_list[-1]) + "\n")

# iterate through data, calculate field values in special cases
for data in clin_data:
    info = data.rstrip('\n').split('\t')
    ovr_surv, age_at_last_known = -1, -1
    sample_to_print = []
    patient_to_print = []
    for i in range(len(info)):
        value = info[i]
        if header[i] == "Last.Known.Clinical.Status":
            if value.startswith("Deceased"):
                value = "DECEASED"
            else:
                value = "LIVING"
        elif header[i] == "Overall.Survival" and value != "Not Reported":
            ovr_surv = int(value)
            value = str(math.floor(float(value)/30.5))
        elif header[i] == "Age.at.Last.Known.Clinical.Status":
            age_at_last_known = int(value)
        elif header[i] == "Progression_Free_Survival":
            if value == "NA":
                if ovr_surv != -1:
                    value = str(math.floor(float(ovr_surv)/30.5))
            else:
                value = str(math.floor(float(value)/30.5))
        if h_dict[header[i]][1] == '1':
            sample_to_print.append(value)
        if h_dict[header[i]][2] == '1':
            patient_to_print.append(value)
    age = float(age_at_last_known)
    if ovr_surv != -1:
        age -= float(ovr_surv) 
    age = str(math.floor(age/365.25))
    patient_to_print.append(age)
    sample_out.write("\t".join(sample_to_print) + "\n")
    patient_out.write("\t".join(patient_to_print) + "\n")
sample_out.close()
patient_out.close()
