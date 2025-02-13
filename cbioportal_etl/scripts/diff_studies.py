#!/usr/bin/env python3
"""
Script to check a CURRENT study on a current cbioportal instance for differences against a local UPDATE build
"""

import argparse
from urllib.parse import urlparse
import glob
import sys
import os
import csv
import typing
import requests


def get_set_venn(left: set[str], right: set[str]) -> tuple[set[str], set[str], set[str]]:
    """Get items in left set not in right + items in right set not in left + items in both sets
    Args:
        left: set of strings
        right: set of strings
    Return:
        left_only: set of strings found in left but not right
        right_only: set of strings found in right but not left
        shared: set of strings shared by left and right
    """
    left_only = left - right
    right_only = right - left
    shared = left & right
    return left_only, right_only, shared

def output_clinical_changes(clin_type: str, current_only_ids: set, update_only_ids: set, current_attr_only: set, update_attr_only: set, attr_cts: dict, suffix: str) -> None:
    """Print change summary to STDOUT and add/delete ids to file based on clin_type
    Args:
        clin_type: str can be PATIENT or SAMPLE
        current_only_ids: set of IDs with info only on current cbioportal instance
        update_only_ids: set of IDs with info only in update flat file
        current_attr_only: set of attributes (basically columns) found only on current cbioportal instance
        update_attr_only: set of attributes (basically columns) found only in update flat file
        attr_cts: dict with counts of attribute changes for clin_type
        suffix: str for output file name suffix
    """
    print( f"{clin_type} CHANGE SUMMARY:" )
    if len(current_only_ids) > 0:
        print( f"{len(current_only_ids)} {clin_type}s in current would be removed" )
        with open(f"delete_id_list_{suffix}", 'w') as del_list:
            writer = csv.writer(del_list, delimiter='\n')
            writer.writerow(current_only_ids)
    if len(update_only_ids) > 0:
        print( f"{len(update_only_ids)} {clin_type}s in update would be added to the current" )
        with open(f"added_id_list_{suffix}", 'w') as add_list:
            writer = csv.writer(add_list, delimiter='\n')
            writer.writerow(update_only_ids)
    if len(current_attr_only) > 0:
        print( f"{len(current_attr_only)} attributes in current would be removed: {','.join(current_attr_only)}" )
    if len(update_attr_only) > 0:
        print( f"{len(update_attr_only)} attributes in update would be added to the current: {','.join(update_attr_only)}" )
    for attr, count in attr_cts.items():
        print( f"{attr} has {count} change(s)" )
    # Print extra newline for readability
    print ("")


def catalog_and_print_delta_ids(current_data: dict, update_data: dict, shared_attrs: set[str], shared_ids: set[str], clin_type: str, out: typing.IO) -> tuple[set, set, set, set, set, dict]:
    """Compare the entries where the current and updated dataset share clinical_ids and attributes. Catalogue the differences (aka deltas) and print them to a file
    Args:
        current_data: dict with data pulled from a current cbioportal instance with clinical id (based on clin_type) as the primary key, clinical attribute as subkey. For example: current_data['7316-7113']['BROAD_HISTOLOGY']='Tumors of sellar region'
        update_data: dict with data pulled from flat file source. Exact format as current_data
        shared_attr: set of attributes shared by current and update
        shared_ids: set of clinical ids shared by current and update
        clin_type: str can be PATIENT or SAMPLE
        out: typing.IO file to print to
    Returns:
        update_delta_ids: set of IDs that exist in both current_data and update_data but have at least 1 differing attribute
        update_delta_attr_cts: dict with counts of attribute changes for clin_type
    """
    # focus on common sample and common attr, as "everything is different for x" is not that useful
    print( f"{clin_type}\tattribute\tbefore\tafter", file=out )
    update_delta_attr_cts: dict = {}
    # track at row level what as changed
    update_delta_ids: set = set()
    for clinical_id in shared_ids:
        for attr in shared_attrs:
            # current will not have a value for that attr in the struct if none
            current_value: str = current_data[clinical_id].get(attr, "NA")
            if current_value != update_data[clinical_id][attr]:
                print( f"{clinical_id}\t{attr}\t{current_value}\t{update_data[clinical_id][attr]}", file=out )
                if attr not in update_delta_attr_cts:
                    update_delta_attr_cts[attr] = 0
                update_delta_attr_cts[attr] += 1
                update_delta_ids.add(clinical_id)
    return update_delta_ids, update_delta_attr_cts


def table_to_dict(in_file: str, key: str, aggr_list: list) -> tuple[dict, set, str, list, str]:
    """Take a text file and convert to dict with certain row value as primary, all other row values as subkeys
    Args:
        in_file: str of input flat file name to parse
        key: str of column name containing primary key to use to aggregate data
        aggr_list: list of column names in which data values are are to be split and sorted to avoid a;b != b;a scenarios
    Returns:
        data_dict: dict with flat file data converted as key value as primary key, associated column name as sub key, and row value as value
        attributes: set of col names sans primary key
        "\t".join(header) +"\n": str col names
        data_clinical: list in_file read in as array of strings
        big_head: str of header lines from in_file

    """
    with open(in_file) as f:
        # track full 5 line header, for diff file printing, uses first non-hash line for finding positions of attributes
        big_head: str = ""
        for entry in f:
            big_head += entry
            if not entry.startswith("#"):
                header: list = entry.rstrip('\n').split('\t')
                break
        data_dict: dict = {}
        data_clinical: list = f.readlines()
        for entry in data_clinical:
            data: list = entry.rstrip('\n').split('\t')
            # Create a dictionary by zipping together the line with the header; sub in "NA" for empty columns
            entry_dict: dict[str, str] = { k: ("NA" if v == "" else v) for k,v in zip(header, data) }
            # Aggregate attributes have multiple entries in a random order separated by a semicolon
            # Therefore, sort them so that when compared, no errors are triggered
            for aggr in aggr_list:
                if aggr in header: entry_dict[aggr] = ';'.join(sorted(entry_dict[aggr].split(';')))
            key_value = entry_dict.pop(key)
            data_dict[key_value] = entry_dict
    attributes: set = set(header)
    # no need for primary key to be reported as an attribute
    attributes.remove(key)
    return data_dict, attributes, "\t".join(header) +"\n", data_clinical, big_head


def data_clinical_from_study(url: str, auth_headers: dict, study_id: str, data_type: str, aggr_list: list) -> dict:
    """Get all the column-value pairs for each data_type(SAMPLE or PATIENT) for a specific study
    Args:
        url: current cbioportal instance url
        auth_headers: dict with authorization token
        study_id: str of study name to pull
        data_type: str PATIENT or SAMPLE
        aggr_list: list of attribute names in which data values are are to be split and sorted to avoid a;b != b;a scenarios
    Returns:
        data_dict: dict with current cbioportal instance data converted as key value as primary key, associated column name as sub key, and row value as value
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
    """Pull all clinical event data (treatment, imaging, etc) for a study currently on a cBio server if the study has such data available
    This will be used to compare to proposed data to update on to the server
    Args:
        url: current cbioportal instance url
        auth_headers: dict with authorization token
        study_id: str of study name to pull
    Returns:
        current_timeline_data: dict of current cbioportal instance timeline data with event type as key
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
    study_timeline_data: list = requests.get(f"{url}/api/studies/{study_id}/clinical-events?projection=DETAILED", headers=auth_headers, timeout=360).json()
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


def data_clinical_timeline_file_to_list(file_path: str) -> tuple[str, list]:
    """
    Read in update flat file headers and data into list
    """
    with open(file_path) as flat_file:
        head = next(flat_file)
        data: list = flat_file.readlines()
    return head, data


def output_delta_by_id(diff_id: set, update_list: list, header: str, id_field: str, outfile_name: str, big_head: str | None) -> None:
    """Use IDs from diff_list, then get all values from update set to output to incremental update data file
    diff_id: set of IDs to filter down update_list for output
    update_list: list of entries from update flat file
    header: str of column names to output unless big_head is populated
    outfile_name: str output file name
    big_head: str lead header lines, if applicable
    """
    id_index: list = header.split('\t').index(id_field)
    update_list = [ item for item in update_list if item.split('\t')[id_index] in diff_id ]
    with open(outfile_name, 'w') as event_delta:
        event_delta.write(header if big_head is None else big_head)
        event_delta.write("".join(update_list))


def compare_timeline_data(current_timeline: dict, update_timeline: dict, out_dir: str, header_info: dict, file_ext: dict) -> None:
    """Compare each event type between current and update candidate to see if any differences exists
    Args:
        current_timeline: dict of current cbioportal instance timeline data with event type as key
        update_timeline: dict of update flat file timeline data with event type as key
        out_dir: sre of output dir location
        header_info: dict of headers by event type for output
        file_ext: dict of file extensions to use by event type for output
    """
    event_type = current_timeline.keys()
    event_ext_dict = {v: k for k, v in file_ext.items()}
    for event in event_type:
        current_ids: set = set(current_timeline[event])
        # temp workaround for testing
        if event not in update_timeline:
            update_timeline[event] = []
        # strip new line
        update_ids = set([ item.rstrip('\n') for item in update_timeline[event] ])
        diff_ids: set = update_ids - current_ids
        if len(diff_ids) > 1:
            print(f"{len(diff_ids)} changes in {event} found for this study. Outputting delta files")
            suffix: str = event_ext_dict[event]
            outfilename: str = f"{out_dir}/data_clinical_timeline_{suffix}"
            uniq_patients: set = [item.split('\t')[0] for item in diff_ids]
            output_delta_by_id(diff_id=uniq_patients, update_list=update_timeline[event], header=header_info[event], id_field="PATIENT_ID", outfile_name=outfilename, big_head=None)


def run_py(args):
    # Create incremental study updates dir
    delta_dir: str = f"{args.study}_delta_data"
    os.makedirs(delta_dir, exist_ok=True)

    # Set Up URL pulling
    print("Processing data clinical info", file=sys.stderr)
    with open(args.token, 'r') as token_file:
        token: str = token_file.read().rstrip().split(': ')[1]
    auth_headers: dict = { "Authorization": f"Bearer {token}"}
    url_object = urlparse(args.url)
    formatted_url: str = f"{url_object.scheme}://{url_object.hostname}"

    # Fetch attribute keys
    fetch_attr_path: str = 'api/clinical-attributes/fetch'
    attr_keys: list[dict[str, str]] = requests.post(f"{formatted_url}/{fetch_attr_path}", headers=auth_headers, json=[args.study], params={"projection":"ID"}, timeout=30).json()

    # Hardcoding a JSON-like payload that will dictate do similar comparisons
    # For each comparison it has implicit and skip attributes
    # Also includes the name of the data table that has the relevant information as well as the key in that table to locate the information
    aggregate_vals: list[str] = ["SPECIMEN_ID", "EXPERIMENT_STRATEGY"]
    comparisons: dict[str, str] = {
        "SAMPLE": {
            "attr_implicit": ['PATIENT_ID'],
            "attr_skip": ['FRACTION_GENOME_ALTERED', 'MUTATION_COUNT'],
            "data_table": args.datasheet_sample,
            "data_table_key": "SAMPLE_ID"
        },
        "PATIENT": {
            "attr_implicit": [],
            "attr_skip": ['SAMPLE_COUNT'],
            "data_table": args.datasheet_patient,
            "data_table_key": "PATIENT_ID"
        }
    }

    for comp in comparisons:
        if comp == "PATIENT":
            current_attr_keys: set[str] = set([x['clinicalAttributeId'] for x in attr_keys if x['patientAttribute']])
        else:
            current_attr_keys: set[str] = set([x['clinicalAttributeId'] for x in attr_keys if not x['patientAttribute']])

        # Remove skip attrs and add implicit attrs
        current_attr_keys -= set(comparisons[comp]["attr_skip"])
        current_attr_keys |= set(comparisons[comp]["attr_implicit"])

        # Pull current data from portal URL; build update data from file
        current_data: dict = data_clinical_from_study(url=formatted_url, auth_headers=auth_headers, study_id=args.study, data_type=comp, aggr_list=aggregate_vals)
        update_data, update_attr_keys, file_header, file_as_list, big_head = table_to_dict(in_file=comparisons[comp]["data_table"], key=comparisons[comp]["data_table_key"], aggr_list=aggregate_vals)

        # Get attrs/ids uniq to the current and update datasets
        current_only_attr: set[str]; update_only_attr: set[str]; shared_attr: set[str]
        current_only_attr, update_only_attr, shared_attr = get_set_venn(left = current_attr_keys, right = update_attr_keys)
        current_only_ids: set[str]; update_only_ids: set[str]; shared_ids: set[str]
        current_only_ids, update_only_ids, shared_ids= get_set_venn(left = set(current_data.keys()), right = set(update_data.keys()))

        # Get and print ids that are present in both but have had updated attrs
        with open(f"{comp.lower()}_current_v_update.txt", 'w') as diff_out:
            update_delta_ids, update_delta_attr_cts = catalog_and_print_delta_ids(current_data=current_data, update_data=update_data, shared_attrs=shared_attr, shared_ids=shared_ids, clin_type=comp, out=diff_out)

        # Get the complete list of new and updated ids
        updated_only_or_delta_ids: set[str] = update_only_ids | update_delta_ids
        if len(updated_only_or_delta_ids) > 1:
            output_delta_by_id( diff_id=updated_only_or_delta_ids, update_list=file_as_list, header=file_header, id_field=comparisons[comp]["data_table_key"], outfile_name=f"{delta_dir}/data_clinical_{comp.lower()}.txt", big_head=big_head )
        output_clinical_changes(clin_type=comp, current_only_ids=current_only_ids , update_only_ids=update_only_ids, current_attr_only=current_only_attr, update_attr_only=update_only_attr, attr_cts=update_delta_attr_cts, suffix=f"{comp}.txt")

    # timeline-diffs
    if args.datasheet_timeline:
        print("Getting current timeline data, please wait...", file=sys.stderr)
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
        # store file data in memory by event type to match current pull dict created
        timeline_update: dict = {}
        # store file headers to repurpose
        event_file_head: dict = {}
        for path in timeline_update_files:
            suffix: str = path.replace(timeline_update_file_prefix, "")
            event_type: str = timeline_event_dict[suffix]
            event_file_head[event_type], timeline_update[event_type] = data_clinical_timeline_file_to_list(file_path=path)
        print("Comparing current timeline to update", file=sys.stderr)
        compare_timeline_data(current_timeline=timeline_current, update_timeline=timeline_update, out_dir=delta_dir, header_info=event_file_head, file_ext=timeline_event_dict)

    print("Done!", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Compare update clinical and timeline data to current cbioportal instance. Outputs changes summaries, id lists, and delta files")
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
