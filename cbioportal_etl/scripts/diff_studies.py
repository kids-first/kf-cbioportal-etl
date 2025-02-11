#!/usr/bin/env python3
"""
Script to check a study on pedcbioportal for differences against a local build
"""

import argparse
from urllib.parse import urlparse
import glob
import sys
import requests
import os
import csv
import typing


def clinical_diffs(current_data: dict, update_data: dict, current_attr: set, update_attr: set, clin_type: str, out: typing.IO, file_header: str, file_as_list: list, out_inc_dir: str, big_head: str) -> None:
    """
    Compare differences in portal sample data and build, output lists for downstream action
    """
    # gross ID diffs
    current_clinical_ids = set(current_data.keys())
    update_clinical_ids = set(update_data.keys())
    # this list is basically rows to be deleted
    current_only_ids: set = current_clinical_ids - update_clinical_ids
    # this list is basically rows to be added
    update_only_ids: set = update_clinical_ids - current_clinical_ids
    common_clinical_ids: set = current_clinical_ids & update_clinical_ids
    # gross attribute diffs
    current_attr_only: set = current_attr - update_attr
    update_attr_only: set = update_attr - current_attr
    common_attr: set = current_attr & update_attr
    # focus on common sample and common attr, as "everything is different for x" is not that useful
    print(clin_type + "\tattribute\tbefore\tafter", file=out)
    attr_cts: dict = {}
    # track at row level what as changed
    id_w_change: set = set()
    for clinical_id in common_clinical_ids:
        for attr in common_attr:
            # portal will not have a value for that attr in the struct if none
            current_value: str = current_data[clinical_id].get(attr, "NA")
            if current_value != update_data[clinical_id][attr]:
                print("{}\t{}\t{}\t{}".format(clinical_id, attr, current_value, update_data[clinical_id][attr]), file=out)
                if attr not in attr_cts:
                    attr_cts[attr] = 0
                attr_cts[attr] += 1
                id_w_change.add(clinical_id)
    # output diff file
    update_ids: set = update_only_ids.union(id_w_change)
    if len(update_ids) > 1:
        suffix = ("sample.txt" if clin_type == "Sample" else "patient.txt")
        outfilename = f"{out_inc_dir}/data_clinical_{suffix}"
        output_delta_by_id( diff_id=update_ids, load_list=file_as_list, header=file_header, id_field=("SAMPLE_ID" if clin_type=="Sample" else "PATIENT_ID"), outfile_name=outfilename, big_head=big_head )
    # print change summary to STDOUT
    print(clin_type +" CHANGE SUMMARY:")
    if len(current_only_ids) > 0:
        print("{} {}s in portal would be removed".format(len(current_only_ids), clin_type))
        with open(f"delete_id_list_{suffix}", 'w') as del_list:
            writer = csv.writer(del_list, delimiter='\n')
            writer.writerow(current_only_ids)
    if len(update_only_ids) > 0:
        print("{} {}s in build would be added to the portal".format(len(update_only_ids), clin_type))
        with open(f"added_id_list_{suffix}", 'w') as add_list:
            writer = csv.writer(add_list, delimiter='\n')
            writer.writerow(update_only_ids)
    if len(current_attr_only) > 0:
        print("{} attributes in portal would be removed: {}".format(len(current_attr_only), ",".join(current_attr_only)))
    if len(update_attr_only) > 0:
        print("{} attributes in build would be added to the portal: {}".format(len(update_attr_only), ",".join(update_attr_only)))
    for attr in attr_cts:
        print("{} has {} change(s)".format(attr, attr_cts[attr]))
    # Print extra newline for readability
    print ("")


def table_to_dict(in_file: str, key: str, aggr_list: list) -> tuple[dict, set, str, list, str]:
    """
    Take a text file and convert to dict with certain row value as primary, all other row values as subkeys.
    Return a set of attribute keys.
    Lastly return file header and list of rows in the event there are load differences
    """
    with open(in_file) as f:
        # track full 5 line header, for diff file printing, uses first non-hash line for finding positions of attributes
        big_head: str = ""
        for entry in f:
            big_head += entry
            if not entry.startswith("#"):
                header: list = entry.rstrip('\n').split('\t')
                primary: str = header.index(key)
                # get aggregate field indices
                aggr_head: list = []
                for aggr in aggr_list:
                    if aggr in header:
                        aggr_head.append(header.index(aggr))
                break
        data_dict: dict = {}
        data_clinical: list = f.readlines()
        for entry in data_clinical:
            data: list = entry.rstrip('\n').split('\t')
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
    attributes: set = set(header)
    # no need for primary key to be reported as an attribute
    attributes.remove(key)
    return data_dict, attributes, "\t".join(header) +"\n", data_clinical, big_head


def data_clinical_from_study(url: str, auth_headers: dict, study_id: str, data_type: str, aggr_list: list) -> dict:
    """
    Get all the column-value pairs for each data_type(SAMPLE or PATIENT) for a specific study
    Convert result to simplified dict
    """
    params: dict = {
        "studyId": study_id,
        "projection": "DETAILED",
        "clinicalDataType": data_type
    }
    get_current_clinical_data_path: str = f"api/studies/{study_id}/clinical-data"
    data_clinical = requests.get(f"{url}/{get_current_clinical_data_path}", headers=auth_headers, params=params, timeout=60).json()
    data_dict: dict = {}
    # Use sampleId or patientID, clinicalAttributeId (column name) and value
    attr_dict: dict = {"SAMPLE": "sampleId", "PATIENT": "patientId" }
    status: list = ["OS_STATUS"]
    for entry in data_clinical:
        clinical_id: str = entry[attr_dict[data_type]]
        if clinical_id not in data_dict:
            # every entry per sample as sampleId and patientId, patient just patientId. Add keys to match
            data_dict[clinical_id] = {"PATIENT_ID": entry['patientId']}
        value: str = entry['value']
        attr_id: str = entry['clinicalAttributeId']
        # For ease of comparison, aggregate attributes concatenated by a separator may not be in the same order, but order does not matter
        # Therefore, sort them so that when compared, no errors are triggered
        if attr_id in aggr_list:
            value: str = ';'.join(sorted(value.split(';')))
        # "standardize" status field so that 0:LIVING = LIVING and 1:DECEASED = DECEASED
        if attr_id in status:
            value: str = value[2:]
        data_dict[clinical_id][attr_id] = value
    return data_dict


def data_clinical_timeline_from_current(url: str, auth_headers: dict, study_id: str) -> dict:
    """
    Pull all clinical event data (treatment, imaging, etc) for a study currently on a cBio server if the study has such data available.
    This will be used to compare to proposed data to load on to the server
    """
    common_fields: list = ["patientId", "startNumberOfDaysSinceDiagnosis", "endNumberOfDaysSinceDiagnosis", "eventType"]
    current_timeline_attr: dict = {
        "Clinical Status": ["CLINICAL_EVENT_TYPE"],
        "IMAGING": ["DIAGNOSTIC_TYPE", "DIAGNOSTIC_TYPE_DETAILED", "BODY_PART", "FLYWHEEL_URL"],
        "SPECIMEN": ["SAMPLE_ID", "SPECIMEN_SITE", "SPECIMEN_TYPE"],
        "SURGERY": ["SUBTYPE", "EXTENT_OF_TUMOR_RESECTION"],
        "TREATMENT": ["TREATMENT_TYPE", "SUBTYPE", "AGENT"]
    }

    current_timeline_data: dict = {}
    for attr in current_timeline_attr:
        current_timeline_data[attr] = []
    clinical_events_study_path: str = f"api/studies/{study_id}/clinical-events?projection=DETAILED"
    study_timeline_data: list = requests.get(f"{url}/{clinical_events_study_path}", headers=auth_headers, timeout=360).json()
    for entry in study_timeline_data:
        event_type: str = entry['eventType']
        temp_event: list = []
        for common in common_fields:
            if common in entry:
                temp_event.append(entry[common])
            else:
                temp_event.append("")
        # event-specific attributes have funky organization,like 'attributes': [{'key': 'CLINICAL_EVENT_TYPE', 'value': 'Initial CNS Tumor'}], easier to flatten it into a normal dict: 
        flatten_attrs: dict = {x['key']: x['value'] for x in entry['attributes']}
        event_specific: list = current_timeline_attr[event_type]
        for field in event_specific:
            if field in flatten_attrs:
                temp_event.append(flatten_attrs[field])
            else:
                temp_event.append("")
        current_timeline_data[event_type].append("\t".join(map(str, temp_event)))
    return current_timeline_data


def data_clinical_timeline_file_to_list(file_path: str) -> (str, list):
    """
    Read in load flat file headers and data into list
    """
    with open(file_path) as flat_file:
        head = next(flat_file)
        data: list = flat_file.readlines()
    return head, data


def output_delta_by_id(diff_id: set, load_list: list, header: str, id_field: str, outfile_name: str, big_head: str | None) -> None:
    """
    Use IDs from diff_list, then get all values from load set to output to incremental load data file
    """
    id_index: list = header.split('\t').index(id_field)
    load_list = [ item for item in load_list if item.split('\t')[id_index] in diff_id ]
    with open(outfile_name, 'w') as event_delta:
        event_delta.write(header if big_head==None else big_head)
        event_delta.write("".join(load_list))


def compare_timeline_data(current_timeline: dict, update_timeline: dict, out_dir: str, header_info: dict, file_ext_dict: dict) -> None:
    """
    Compare each event type between portal and load candidate to see if any differences exists
    """
    event_type = current_timeline.keys()
    event_ext_dict = {v: k for k, v in file_ext_dict.items()}
    for event in event_type:
        current_ids: set = set(current_timeline[event])
        # strip new line
        update_ids = set([ item.rstrip('\n') for item in update_timeline[event] ])
        diff_ids: set = update_ids - current_ids
        if len(diff_ids) > 1:
            print(f"{len(diff_ids)} changes in {event} found for this study. Outputting delta files")
            suffix: str = event_ext_dict[event]
            outfilename: str = f"{out_dir}/data_clinical_timeline_{suffix}"
            uniq_patients: set = [item.split('\t')[0] for item in diff_ids]
            output_delta_by_id(diff_id=uniq_patients, load_list=update_timeline[event], header=header_info[event], id_field="PATIENT_ID", outfile_name=outfilename, big_head=None)


def run_py(args):
    with open(args.token, 'r') as token_file:
        token: str = token_file.read().rstrip().split(': ')[1]

    url_object = urlparse(args.url)
    formatted_url: str = f"{url_object.scheme}://{url_object.hostname}"
    auth_headers: dict = { "Authorization": f"Bearer {token}"}

    # Create incremental study updates dir
    update_dir: str = f"{args.study}_data"
    os.makedirs(update_dir, exist_ok=True)
    # hardcode for now names of aggregate fields, implicit, and skip fields
    aggregate_vals: list = ["SPECIMEN_ID", "EXPERIMENT_STRATEGY"]
    current_sample_attr_implicit: list = ['PATIENT_ID']
    current_patient_attr_skip: list = ['SAMPLE_COUNT']
    current_sample_attr_skip: list = ['FRACTION_GENOME_ALTERED', 'MUTATION_COUNT']

    # get attribute keys
    print("Processing data clinical info", file=sys.stderr)
    fetch_attr_path: str = 'api/clinical-attributes/fetch'
    attr_keys: list = requests.post(f"{formatted_url}/{fetch_attr_path}", headers=auth_headers, json=[args.study], params={"projection":"ID"}, timeout=30).json()
    current_sample_attr_keys: set = set([x['clinicalAttributeId'] for x in attr_keys if not x['patientAttribute']])
    current_patient_attr_keys: set = set([x['clinicalAttributeId'] for x in attr_keys if x['patientAttribute']])

    # gather sample-level metadata
    current_sample_data: dict = data_clinical_from_study(url=formatted_url, auth_headers=auth_headers, study_id=args.study, data_type="SAMPLE", aggr_list=aggregate_vals)
    update_sample_data, update_sample_attr_keys, sample_file_header, sample_file_as_list, sample_big_head = table_to_dict(in_file=args.datasheet_sample, key="SAMPLE_ID", aggr_list=aggregate_vals)
    # implicit attributes not returned by function that are required for sample view
    current_sample_attr_keys.update(current_sample_attr_implicit)
    # drop attributes that are post-load portal-specific
    current_sample_attr_keys -= set(current_sample_attr_skip)
    # sample-level diffs
    with open('sample_portal_v_build.txt', 'w') as sample_diff_out:
        clinical_diffs(current_data=current_sample_data, update_data=update_sample_data, current_attr=current_sample_attr_keys, update_attr=update_sample_attr_keys, clin_type="Sample", out=sample_diff_out, file_header=sample_file_header, file_as_list=sample_file_as_list, out_inc_dir=update_dir, big_head=sample_big_head)
    # patient-level diffs
    current_patient_data: dict = data_clinical_from_study(url=formatted_url, auth_headers=auth_headers, study_id=args.study, data_type="PATIENT", aggr_list=aggregate_vals)
    update_patient_data, update_patient_attr_keys, patient_file_header, patient_file_as_list, patient_big_head = table_to_dict(in_file=args.datasheet_patient, key="PATIENT_ID", aggr_list=aggregate_vals)
    # timeline-diffs
    if args.datasheet_timeline:
        print("Getting portal timeline data, please wait...", file=sys.stderr)
        timeline_event_dict: dict = {
            "clinical_event.txt": "Clinical Status",
            "imaging.txt": "IMAGING",
            "specimen.txt": "SPECIMEN",
            "surgery.txt": "SURGERY",
            "treatment.txt": "TREATMENT", 
        }
        timeline_current: dict = data_clinical_timeline_from_current(url=formatted_url, auth_headers=auth_headers, study_id=args.study)
        timeline_dir: str = args.datasheet_timeline.rstrip("/")
        timeline_update_file_prefix: str = f"{timeline_dir}/data_clinical_timeline_"
        timeline_update_files: list = glob.glob(f"{timeline_update_file_prefix}*")
        # store file data in memory by event type to match portal pull dict created
        timeline_update: dict = {}
        # store file headers to repurpose
        event_file_head: dict = {}
        for path in timeline_update_files:
            suffix: str = path.replace(timeline_update_file_prefix, "")
            event_type: str = timeline_event_dict[suffix]
            event_file_head[event_type], timeline_update[event_type] = data_clinical_timeline_file_to_list(file_path=path)
        print("Comparing current timeline to update", file=sys.stderr)
        compare_timeline_data(current_timeline=timeline_current, update_timeline=timeline_update, out_dir=update_dir, header_info=event_file_head, file_ext_dict=timeline_event_dict)

    # drop attributes that are post-load portal-specific
    current_patient_attr_keys -= set(current_patient_attr_skip)
    with open('patient_portal_v_build.txt', 'w') as patient_diff_out:
        clinical_diffs(current_data=current_patient_data, update_data=update_patient_data, current_attr=current_patient_attr_keys, update_attr=update_patient_attr_keys, clin_type="Patient", out=patient_diff_out, file_header=patient_file_header, file_as_list=patient_file_as_list, out_inc_dir=update_dir, big_head=patient_big_head)
    print("Done!", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Compare local clinical data to server")
    parser.add_argument("-u", "--url", action="store", dest="url", help="url to search against", default="https://pedcbioportal.kidsfirstdrc.org/api/v2/api-docs")
    parser.add_argument("-s", "--study", action="store", dest="study", help="Cancer study ID to compare on server")
    parser.add_argument("-t", "--token", action="store", dest="token", help="Token file obtained from Web API")
    parser.add_argument("-ds", "--datasheet-sample", action="store", dest="datasheet_sample", help="File containing cBio-formatted sample metadata, typically named data_clinical_sample.txt")
    parser.add_argument("-dp", "--datasheet-patient", action="store", dest="datasheet_patient", help="File containing cBio-formatted patient metadata, typically named data_clinical_patient.txt")
    parser.add_argument("-dt", "--datasheet-timeline", action="store", dest="datasheet_timeline", help="Dir containing cBio-formatted timeline data, typically named data_clinical_timeline*")

    args = parser.parse_args()
    run_py(args)


if __name__ == '__main__':
    main()
