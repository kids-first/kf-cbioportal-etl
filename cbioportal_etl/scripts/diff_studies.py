#!/usr/bin/env python3
"""
Script to check a study on pedcbioportal for differences against a local build
"""

import argparse
from urllib.parse import urlparse
import requests
import glob
import pdb


def clinical_diffs(portal, build, portal_attr, build_attr, clin_type, out):
    """
    Compare differences in portal sample data and build.
    """
    # gross ID diffs
    portal_clinical_ids = set(portal.keys())
    build_clinical_ids = set(build.keys())
    portal_only = sorted(portal_clinical_ids - build_clinical_ids)
    build_only = sorted(build_clinical_ids - portal_clinical_ids)
    common_clinical_ids = sorted(portal_clinical_ids & build_clinical_ids)
    # gross attribute diffs
    portal_attr_only = list(portal_attr - build_attr)
    build_attr_only = list(build_attr - portal_attr)
    common_attr = list(portal_attr & build_attr)
    # focus on common samp and common attr, as "everything is different for x" is not that useful
    print(clin_type + "\tattribute\tbefore\tafter", file=out)
    attr_cts = {}
    for clinical_id in common_clinical_ids:
        for attr in common_attr:
            # portal will not have a value for that attr in the struct if none
            portal_value = portal[clinical_id].get(attr, "NA")
            if portal_value != build[clinical_id][attr]:
                print("{}\t{}\t{}\t{}".format(clinical_id, attr, portal_value, build[clinical_id][attr]), file=out)
                if attr not in attr_cts:
                    attr_cts[attr] = 0
                attr_cts[attr] += 1

    # print change summary to STDOUT
    print(clin_type +" CHANGE SUMMARY:")
    if len(portal_only) > 0:
        print("{} {}s in portal would be removed: {}".format(len(portal_only), clin_type, ",".join(portal_only)))
    if len(build_only) > 0:
        print("{} {}s in build would be added to the portal: {}".format(len(build_only), clin_type,  ",".join(build_only)))
    if len(portal_attr_only) > 0:
        print("{} attributes in portal would be removed: {}".format(len(portal_attr_only), ",".join(portal_attr_only)))
    if len(build_attr_only) > 0:
        print("{} attributes in build would be added to the portal: {}".format(len(build_attr_only), ",".join(build_attr_only)))
    for attr in attr_cts:
        print("{} has {} change(s)".format(attr, attr_cts[attr]))
    # Print extra newline for readability
    print ("")


def table_to_dict(in_file, key, aggr_list):
    """
    Take a text file and convert to dict with certain row value as primary, all other row values as subkeys.
    Also return a set of attribute keys
    """
    with open(in_file) as f:
        # skip lines starting with hash until normal header is reached
        for entry in f:
            if not entry.startswith("#"):
                header = entry.rstrip('\n').split('\t')
                primary = header.index(key)
                # get aggregate field indices
                aggr_head = []
                for aggr in aggr_list:
                    if aggr in header:
                        aggr_head.append(header.index(aggr))
                break
        data_dict = {}
        for entry in f:
            data = entry.rstrip('\n').split('\t')
            # Replace empty string with NA as that is how the portal will return it
            data = ["NA" if d == "" else d for d in data]
            data_dict[data[primary]] = {}
            # For ease of comparison, aggregate attributes concatenated by a separator may not be in the same order, but order does not matter
            # Therefore, sort them so that when compared, no errors are triggered
            for i in aggr_head:
                data[i] = ';'.join(sorted(data[i].split(';')))
            # two loops, for up until primary key, then after.
            for i in range(len(data)):
                if i == primary: continue
                data_dict[data[primary]][header[i]] = data[i]
    attr_set = set(header)
    # no need for primary key to be reported as an attribute
    attr_set.remove(key)
    return data_dict, attr_set


def data_clinical_from_study(url, auth_headers, study_id, data_type, aggr_list):
    """
    Get all the column-value pairs for each data_type(SAMPLE or PATIENT) for a specific study
    Convert result to simplified dict
    """
    params = {
        "studyId": study_id,
        "projection": "DETAILED",
        "clinicalDataType": data_type
    }
    get_all_clinical_data_study_path = f"api/studies/{study_id}/clinical-data"
    data_clinical = requests.get(f"{url}/{get_all_clinical_data_study_path}", headers=auth_headers, params=params, timeout=60).json()
    data_dict = {}
    # Use sampleId or patientID, clinicalAttributeId (column name) and value
    attr_dict = {"SAMPLE": "sampleId", "PATIENT": "patientId" }
    status = ["OS_STATUS"]
    for entry in data_clinical:
        clinical_id = entry[attr_dict[data_type]]
        if clinical_id not in data_dict:
            # every entry per sample as sampleId and patientId, patient just patientId. Add keys to match
            data_dict[clinical_id] = {"PATIENT_ID":entry['patientId']}
        value = entry['value']
        attr_id = entry['clinicalAttributeId']
        # For ease of comparison, aggregate attributes concatenated by a separator may not be in the same order, but order does not matter
        # Therefore, sort them so that when compared, no errors are triggered
        if attr_id in aggr_list:
            value = ';'.join(sorted(value.split(';')))
        # "standardize" status field so that 0:LIVING = LIVING and 1:DECEASED = DECEASED
        if attr_id in status:
            value = value[2:]
        data_dict[clinical_id][attr_id] = value
    return data_dict


def data_clinical_treatment_from_study(url, auth_headers, study_id):
    """
    Get and compare treatment data when available
    """
    # from https://realpython.com/python-data-structures/#dataclassesdataclass-python-37-data-classes
    common_fields = ["patientId", "startNumberOfDaysSinceDiagnosis", "endNumberOfDaysSinceDiagnosis", "eventType"]
    portal_timeline_attr_dict = {
        "Clinical Status": ["CLINICAL_EVENT_TYPE"],
        "IMAGING": ["DIAGNOSTIC_TYPE", "DIAGNOSTIC_TYPE_DETAILED", "BODY_PART", "FLYWHEEL_URL"],
        "SPECIMEN": ["SAMPLE_ID", "SPECIMEN_SITE", "SPECIMEN_TYPE"],
        "SURGERY": ["EVENT_TYPE", "SUBTYPE", "EXTENT_OF_TUMOR_RESECTION"],
        "TREATMENT": ["TREATMENT_TYPE", "SUBTYPE", "AGENT"]
    }

    clinical_events_study_path = f"api/studies/{study_id}/clinical-events?projection=DETAILED"
    study_timeline_data = requests.get(f"{url}/{clinical_events_study_path}", headers=auth_headers, timeout=360).json()
    portal_timeline_data_dict = {}
    for attr in portal_timeline_attr_dict:
        portal_timeline_data_dict[attr] = []
    for entry in study_timeline_data:
        event_type = entry['eventType']
        temp_event_list = []
        for common in common_fields:
            if common in entry:
                temp_event_list.append(entry[common])
            else:
                temp_event_list.append("")
        event_specific = portal_timeline_attr_dict[event_type]
        # event-specific attributes have funky organization,like 'attributes': [{'key': 'CLINICAL_EVENT_TYPE', 'value': 'Initial CNS Tumor'}], easier to flatten it into a normal dict: 
        flatten_attrs = {k: v for d in [ {x['key']: x['value']} for x in entry['attributes'] ] for k, v in d.items()}
        for field in event_specific:
            if field in flatten_attrs:
                temp_event_list.append(flatten_attrs[field])
            else:
                temp_event_list.append("")
        portal_timeline_data_dict[event_type].append("\t".join(map(str, temp_event_list)))


    return portal_timeline_data_dict


def data_clinical_timeline_file_to_dict(file_path):
    with open(file_path) as flat_file:
        data_list = flat_file.read().splitlines()
    return data_list


def compare_timeline_data(portal_timeline, load_timeline):
    event_type = portal_timeline.keys()
    for event in event_type:
        portal_set = set(portal_timeline[event])
        load_set = set(load_timeline[event])
        print(load_set - portal_set)


def run_py(args):
    with open(args.token, 'r') as token_file:
        token = token_file.read().rstrip().split(': ')[1]

    url_object = urlparse(args.url)
    formatted_url = "{}://{}".format(url_object.scheme, url_object.hostname)
    auth_headers = { "Authorization": f"Bearer {token}"}
 
    # hardcode for now names of aggregate fields, implicit, and skip fields
    aggr_list = ["SPECIMEN_ID", "EXPERIMENT_STRATEGY"]
    portal_sample_attr_implicit = ['PATIENT_ID']
    portal_patient_attr_skip = ['SAMPLE_COUNT']
    portal_sample_attr_skip = ['FRACTION_GENOME_ALTERED', 'MUTATION_COUNT']

    # get attribute keys
    fetch_attr_path = 'api/clinical-attributes/fetch'
    attr_key_json = requests.post(f"{formatted_url}/{fetch_attr_path}", headers=auth_headers, json=[args.study], params={"projection":"ID"}, timeout=30).json()
    portal_sample_attr_keys = set([x['clinicalAttributeId'] for x in attr_key_json if not x['patientAttribute']])
    portal_patient_attr_keys = set([x['clinicalAttributeId'] for x in attr_key_json if x['patientAttribute']])

    # gather sample-level metadata
    portal_sample_data = data_clinical_from_study(url=formatted_url, auth_headers=auth_headers, study_id=args.study, data_type="SAMPLE", aggr_list=aggr_list)
    build_sample_data, build_sample_attr_keys = table_to_dict(args.datasheet_sample, "SAMPLE_ID", aggr_list)
    # implicit attributes not returned by function that are required for sample view
    portal_sample_attr_keys.update(portal_sample_attr_implicit)
    # drop attributes that are post-load portal-specific
    portal_sample_attr_keys -= set(portal_sample_attr_skip)
    # sample-level diffs
    with open('sample_portal_v_build.txt', 'w') as sample_diff_out:
        clinical_diffs(portal_sample_data, build_sample_data, portal_sample_attr_keys, build_sample_attr_keys, "Sample", sample_diff_out)
    # patient-level diffs
    portal_patient_data = data_clinical_from_study(url=formatted_url, auth_headers=auth_headers, study_id=args.study, data_type="PATIENT", aggr_list=aggr_list)
    build_patient_data, build_patient_attr_keys = table_to_dict(args.datasheet_patient, "PATIENT_ID", aggr_list)
    # timeline-diffs
    if args.datasheet_timeline:
        timeline_event_dict = {
            "clinical_event.txt": "Clinical Status",
            "imaging.txt": "IMAGING",
            "specimen.txt": "SPECIMEN",
            "surgery.txt": "SURGERY",
            "treatment.txt": "TREATMENT", 
        }
        timeline_portal = data_clinical_treatment_from_study(url=formatted_url, auth_headers=auth_headers, study_id=args.study)
        timeline_dir = args.datasheet_timeline.rstrip("/")
        timeline_flatfile_prefix = f"{timeline_dir}/data_clinical_timeline_"
        timeline_flatfile_list = glob.glob(f"{timeline_flatfile_prefix}*")
        timeline_flatfile_dict = {}
        for path in timeline_flatfile_list:
            suffix = path.replace(timeline_flatfile_prefix, "")
            timeline_flatfile_dict[timeline_event_dict[suffix]] = data_clinical_timeline_file_to_dict(path)
        compare_timeline_data(portal_timeline=timeline_portal, load_timeline=timeline_flatfile_dict)

    # drop attributes that are post-load portal-specific
    portal_patient_attr_keys -= set(portal_patient_attr_skip)
    with open('patient_portal_v_build.txt', 'w') as patient_diff_out:
        clinical_diffs(portal_patient_data, build_patient_data, portal_patient_attr_keys, build_patient_attr_keys, "Patient", patient_diff_out)


def main():
    parser = argparse.ArgumentParser(description="Compare local clinical data to server")
    parser.add_argument("-u", "--url", action="store", dest="url", help="url to search against", default="https://pedcbioportal.kidsfirstdrc.org")
    parser.add_argument("-c", "--study", action="store", dest="study", help="Cancer study ID to compare on server")
    parser.add_argument("-t", "--token", action="store", dest="token", help="Token file obtained from Web API")
    parser.add_argument("-ds", "--datasheet-sample", action="store", dest="datasheet_sample", help="File containing cBio-formatted sample metadata, typically named data_clinical_sample.txt")
    parser.add_argument("-dp", "--datasheet-patient", action="store", dest="datasheet_patient", help="File containing cBio-formatted patient metadata, typically named data_clinical_patient.txt")
    parser.add_argument("-dt", "--datasheet-timeline", action="store", dest="datasheet_timeline", help="Dir containing cBio-formatted timeline data, typically named data_clinical_timeline*")

    args = parser.parse_args()
    run_py(args)


if __name__ == '__main__':
    main()
