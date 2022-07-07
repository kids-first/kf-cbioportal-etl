import sys
import argparse
import json
import math
import re
import pdb


def format_desc(entry):
    desc = entry.replace("_", " ")
    return desc

def format_colname(entry):
    col_name = entry.replace(".", "_").upper()
    return col_name


def init_list(header_list, entry):
    '''
    Initialize list of list of header entries
    '''
    header_list.append([])
    if h_dict[entry][m_idx] == "":
        header_list[0].append(entry)
    else:
        header_list[0].append(h_dict[entry][m_idx])
    header_list.append([])
    if h_dict[entry][d_idx] == "":
        header_list[1].append(format_desc(entry))
    else:
        header_list[1].append(h_dict[entry][d_idx])
    header_list.append([])
    header_list[2].append(h_dict[entry][t_idx])
    header_list.append([])
    header_list[3].append(h_dict[entry][o_idx])
    header_list.append([])
    if h_dict[entry][0] == "":
        header_list[4].append(format_colname(entry))
    else:
        header_list[4].append(h_dict[entry][0])
    return header_list


def build_header(header_list, entry):
    '''
    Keep adding to existing entries in header list of list
    '''
    if not header_list:
        header_list = init_list(header_list, entry)
    else:
        if h_dict[entry][m_idx] == "":
            header_list[0].append(entry)
        else:
            header_list[0].append(h_dict[entry][m_idx])
        if h_dict[entry][d_idx] == "":
            header_list[1].append(format_desc(entry))
        else:
            header_list[1].append(h_dict[entry][d_idx])
        header_list[2].append(h_dict[entry][t_idx])
        header_list[3].append(h_dict[entry][o_idx])
        if h_dict[entry][m_idx] == "":
            header_list[4].append(format_colname(entry))
        else:
            header_list[4].append(h_dict[entry][m_idx])

    return header_list

parser = argparse.ArgumentParser(description='Script to convert clinical data to cbio clinical data sheets')
parser.add_argument('-f', '--header-file', action='store', dest="head" ,help='tsv file with input file original sample '
                    'names, output sheet flag, and conversion')
parser.add_argument('-c', '--clinical-data', action='store', dest='clin',
                    help='Input clinical data sheet')
parser.add_argument('-s', '--cell-line-supplement', action='store', dest='cl_supp', help='supplemental file with cell '
                                                                        'line meta data - bs id<tab>type.  optional')


args = parser.parse_args()

# read in header key file to assign some of the clinical data headers to cbio standard ones and to all a data type and sheet location
# File header is: original	MAP TO	sample	patient	type
h_dict = {}
h_file = open(args.head)
h_head = next(h_file)
h_header = h_head.rstrip('\n').split('\t')
# Get key header loc, shift-1 since fist col used as key
# header cols currently: original	MAP TO	description	sample	patient	type	order
m_idx = h_header.index("MAP TO") - 1
d_idx = h_header.index("description") - 1
s_idx = h_header.index("sample") - 1
p_idx = h_header.index("patient") - 1
t_idx = h_header.index("type") - 1
o_idx = h_header.index("order") - 1
for line in h_file:
    info = line.rstrip('\n').split('\t')
    h_dict[info[0]] = info[1:]
# store cell line data
cl_supp = {}
if args.cl_supp is not None:
    cl_info = open(args.cl_supp)
    for line in cl_info:
        info = line.rstrip('\n').split('\t')
        cl_supp[info[0]] = info[1]


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
sample_out = open("data_clinical_sample.txt", "w")
patient_out = open("data_clinical_patient.txt", "w")
# initialize first entry, as it has all leading "#""
patient_head_list = []
sample_head_list = []
for entry in header:
    if entry in h_dict:
        if h_dict[entry][s_idx] == '1':
            sample_head_list = build_header(sample_head_list, entry)
        if h_dict[entry][p_idx] == '1':
            patient_head_list = build_header(patient_head_list, entry)


# output headers to files
for i in range(0, 4, 1):
    sample_out.write("#" + "\t".join(sample_head_list[i]) + "\n")
    patient_out.write("#" + "\t".join(patient_head_list[i]) + "\n")
sample_out.write("\t".join(sample_head_list[-1]) + "\n")
patient_out.write("\t".join(patient_head_list[-1]) + "\n")
# pt dict to track if already seen to avoid duplicate output
pt_id_dict = {}
# convert values to those that have function and meaning in cbio
tumor_descriptor_dict = {
    "Initial CNS Tumor": "primary", 
    "Primary Tumor": "primary", "Metastatic Tumor": "metastasis",
    "Second Malignancy": "second malignancy", "Recurrence": "recurrence",
    "Progressive": "progressed",
    "Progressive Disease Post Mortem": "progressed"
    }

# track some sample info so that somatic events can be collapsed
samp_dict = {}
s_type = header.index('sample_type')
comp = header.index('composition')
bs_id = header.index('Kids_First_Biospecimen_ID')
exp = header.index('experimental_strategy')
id_mapping = {}
bs_type = {}
# iterate through data, calculate field values in special cases

for data in clin_data:
    info = data.rstrip('\n').split('\t')
    if info[s_type] == "Normal" and info[exp] != 'RNA-Seq':
        continue
    pt_id = info[header.index("Kids_First_Participant_ID")]
    if pt_id in pt_id_dict:
        pt_id_dict[pt_id] += 1
    else:
        pt_id_dict[pt_id] = 1

    sample_to_print = []
    patient_to_print = []
    samp_id = ""
    for i in range(len(info)):
        if header[i] in h_dict:
            value = info[i]
            if h_dict[header[i]][t_idx] == "NUMBER":
                try:
                    float(value)
                except ValueError:  
                    value = "NA"
                if header[i] == "OS_days" and value != "NA" and value != "":
                    value = str(math.floor(float(value)/30.5))
                    # ovr_surv = value
                elif header[i] == "age_at_diagnosis_days":
                    # age_at_last_known = value
                    if value != "NA":
                        # age_at_last_known = float(value)
                        # age_at_last_known = str(math.floor(age_at_last_known/365.25))
                        value = str(math.floor(float(value)/365.25))
                elif header[i] == "PFS_days" and value != "NA":
                    value = str(math.floor(float(value)/30.5))
                    # d_free_mos = value
            elif header[i] == "tumor_descriptor":
                if value in tumor_descriptor_dict:
                    value = tumor_descriptor_dict[value]
            # replace status with NA if value not acceptable
            elif header[i] == "OS_status":
                if value not in ["LIVING", "DECEASED", "NA"]:
                    sys.stderr.write("WARN: OS_status was " + value + ", setting to NA\n")
                    value = "NA"
            # Keep track of sample id for QC, etc
            elif header[i] == "formatted_sample_id":
                samp_id = value
            if h_dict[header[i]][s_idx] == '1':
                sample_to_print.append(value)
            if h_dict[header[i]][p_idx] == '1' and pt_id_dict[pt_id] == 1:
                patient_to_print.append(value)
    if pt_id_dict[pt_id] == 1:
        patient_out.write("\t".join(patient_to_print) + "\n")
    # sample_out.write("\t".join(sample_to_print) + "\n")
    if samp_id not in samp_dict:
        samp_dict[samp_id] = sample_to_print
        id_mapping[samp_id] = []
    id_mapping[samp_id].append(info[bs_id])
    if info[exp] == "RNA-Seq":
        bs_type[info[bs_id]] = "RNA"
    else:
        bs_type[info[bs_id]] = "DNA"
# cycle through sample IDs to see if there's matched DNA/RNA and if one can be made
check = {}
for samp_id in id_mapping:
    if len(id_mapping[samp_id]) == 2:
        # QC check, should only be one DNA and one RNA
        if bs_type[id_mapping[samp_id][0]] == bs_type[id_mapping[samp_id][1]]:
            sys.stderr.write("Duplicate assay types for  " + samp_id + ": " + ",".join(id_mapping[samp_id]) + "\n")
            # exit(1)
        spec = id_mapping[samp_id][0] + ";" + id_mapping[samp_id][1]
        if bs_type[id_mapping[samp_id][0]] == "RNA":
            spec = id_mapping[samp_id][1] + ";" + id_mapping[samp_id][0]
        samp_dict[samp_id][0] = spec
    elif len(id_mapping[samp_id]) > 2:
        # QC check, only one or two biospec per sample ID
        sys.stderr.write("Saw more than two biospecimens for " + samp_id + ": " + ",".join(id_mapping[samp_id]) + "\n")
        # exit(1)
    else:
        # skip cell line re-matching
        if not re.search("-CL", samp_id):
            parts = samp_id.split("_")
            if parts[0] not in check:
                check[parts[0]] = []
            check[parts[0]].append(samp_id)

for samp_id in samp_dict:
    sample_out.write("\t".join(samp_dict[samp_id]) + "\n")
mapping_file = open("bs_id_sample_map.txt", "w")
mapping_file.write("BS_ID\tSample Type\tCbio ID\n")
for samp_id in id_mapping:
    for bs_id in id_mapping[samp_id]:
        try:
            mapping_file.write("\t".join([bs_id, bs_type[bs_id], samp_id]) + "\n")
        except Exception as e:
            sys.stderr.write(str(e) + "\n")
            pdb.set_trace()
            hold = 1
sample_out.close()
patient_out.close()
mapping_file.close()
