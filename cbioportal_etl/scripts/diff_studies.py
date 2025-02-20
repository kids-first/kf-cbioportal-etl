#!/usr/bin/env python3
"""Script to check a CURRENT study on a current cbioportal instance for differences against a local UPDATE build.

Outputs files that are the foundation for incremental updates and clinical data QC
"""

import argparse
import csv
import glob
import os
import sys
import typing
from urllib.parse import urlparse
from typing import TypedDict

import requests


def parse_file(
    file_path: str, header_symbol: str | None = None, column_names: bool = False
) -> tuple[list[str] | None, list[str] | None, list[str]]:
    """Break up a file into header, column names, and body

    Args:
        file_path: string path to file to process
        column_names: boolean value declaring whether the file has a line that contains the column names (first line after header)
        header_symbol: string value that precedes all
    Return:
        header: List of strings containing all lines with a header symbol and the column name line
        columns: List of strings containing the column name line split by tab
        body: List of strings containing all non-header lines in the file
    """
    with open(file_path, "rt") as fh:
        all_lines = fh.readlines()
    end_header_index = 0
    if header_symbol:
        for index, line in enumerate(all_lines):
            if not line.startswith(header_symbol):
                end_header_index = index
                break
    columns = None
    if column_names:
        columns = all_lines[end_header_index].strip().split("\t")
        end_header_index += 1
    header = all_lines[:end_header_index] if len(all_lines[:end_header_index]) > 0 else None
    body = all_lines[end_header_index:]
    return header, columns, body


def print_parsed_file(header: list[str] | None, body: list[str], out_filename: str) -> None:
    with open(out_filename, "wt") as out_file:
        if header:
            out_file.write("".join(header))
        out_file.write("".join(body))


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


def output_clinical_changes(
    clin_type: str,
    current_only_ids: set,
    update_only_ids: set,
    current_attr_only: set,
    update_attr_only: set,
    attr_cts: dict,
    suffix: str,
) -> None:
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
    print(f"{clin_type} CHANGE SUMMARY:")
    if len(current_only_ids) > 0:
        print(f"{len(current_only_ids)} {clin_type}s in current would be removed")
        with open(f"delete_id_list_{suffix}", "w") as del_list:
            writer = csv.writer(del_list)
            writer.writerow(current_only_ids)
    if len(update_only_ids) > 0:
        print(f"{len(update_only_ids)} {clin_type}s in update would be added to the current")
    if len(current_attr_only) > 0:
        print(
            f"{len(current_attr_only)} attributes in current would be removed: {','.join(current_attr_only)}"
        )
    if len(update_attr_only) > 0:
        print(
            f"{len(update_attr_only)} attributes in update would be added to the current: {','.join(update_attr_only)}"
        )
    for attr, count in attr_cts.items():
        print(f"{attr} has {count} change(s)")
    # Print extra newline for readability
    print("")


def split_data_file(
        data_lines: list[str],
        column_names: list[str],
        key: str,
        aggr_list: list,
        current_data: dict[str,str],
        shared_attrs: set[str],
) -> tuple[dict[str, dict[str, str]], list[str], dict[str, dict[str,str]], list[str] | None]:
    """Take a text file and convert to dict with certain row value as primary, all other row values as subkeys

    Args:
        data_lines: list of strings representing the lines of data
        column_names: list of strings representing the column names
        key: str of column name to use as primary key in dict for all other columns
        aggr_list: list of column names in which data values are to be split and sorted to avoid a;b != b;a scenarios
        current_data: dict containing the current dataset
        shared_attrs: set containing the current attributes
    Returns:
        delta_data
        delta_lines
        new_data
        new_lines
    """
    delta_data: dict[str, dict] = {}
    delta_lines: list[str] = []
    new_data: dict[str, dict] = {}
    new_lines: list[str] = []
    for entry in data_lines:
        data: list = entry.rstrip("\n").split("\t")
        # Create a dictionary by zipping together the line with the header; sub in "NA" for empty columns
        entry_dict: dict[str, str] = {k: ("NA" if v == "" else v) for k, v in zip(column_names, data)}
        # Aggregate attributes have multiple entries in a random order separated by a semicolon
        # Therefore, sort them so that when compared, no errors are triggered
        for aggr in aggr_list:
            if aggr in header:
                entry_dict[aggr] = ";".join(sorted(entry_dict[aggr].split(";")))
        key_value = entry_dict.pop(key)
        if key_value in current_data and is_delta(current_data[key_value], entry_dict, shared_attrs):
            delta_data[key_value] = entry_dict
            delta_lines.append(entry)
        else:
            new_data[key_value] = entry_dict
            new_lines.append(entry)
    return delta_data, delta_lines, new_data, new_lines


def is_delta(current_data, update_data, shared_attrs) -> bool:
    """Given two dicts, compare their shared attributes and report if any of their attributes differ

    :param current_data:
    :param update_data:
    :param shared_attrs:
    :return:
    """
    for attr in shared_attrs:
        # current will not have a value for that attr in the struct if none
        current_value: str = current_data.get(attr, "NA")
        if current_value != update_data[attr]:
            return True
    return False

def print_and_count_delta_attr(
        clin_type: str,
        delta_data: dict[str,str],
        current_data: dict[str, str],
        shared_attr: set[str],
        out_filename: str
) -> dict[str, int]:
    delta_attr_cts: dict[str, int] = {}
    with open(out_filename, "w") as diff_out:
        print(f"{clin_type}\tattribute\tbefore\tafter", file=diff_out)
        for clinical_id in delta_data:
            for attr in shared_attr:
                current_value: str = current_data[clinical_id].get(attr, "NA")
                if current_value != delta_data[clinical_id][attr]:
                    print(
                        f"{clinical_id}\t{attr}\t{current_value}\t{delta_data[clinical_id][attr]}",
                        file=diff_out,
                    )
                    if attr not in update_delta_attr_cts:
                        delta_attr_cts[attr] = 0
                    delta_attr_cts[attr] += 1
    return delta_attr_cts

def data_clinical_from_study(
    url: str, auth_headers: dict, study_id: str, data_type: str, aggr_list: list
) -> dict:
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
        "clinicalDataType": data_type,
    }
    get_current_clinical_data_path: str = f"api/studies/{study_id}/clinical-data"
    data_clinical = requests.get(
        f"{url}/{get_current_clinical_data_path}",
        headers=auth_headers,
        params=params,
        timeout=60,
    ).json()
    data_dict: dict = {}
    # Use sampleId or patientID, clinicalAttributeId (column name) and value
    attr_dict: dict = {"SAMPLE": "sampleId", "PATIENT": "patientId"}
    status: list = ["OS_STATUS"]
    for entry in data_clinical:
        clinical_id: str = entry[attr_dict[data_type]]
        if clinical_id not in data_dict:
            # every entry per sample as sampleId and patientId, patient just patientId. Add keys to match
            data_dict[clinical_id] = {"PATIENT_ID": entry["patientId"]}
        value: str = entry["value"]
        attr_id: str = entry["clinicalAttributeId"]
        # For ease of comparison, aggregate attributes concatenated by a separator may not be in the same order, but order does not matter
        # Therefore, sort them so that when compared, no errors are triggered
        if attr_id in aggr_list:
            value: str = ";".join(sorted(value.split(";")))
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
    common_fields: list = [
        "patientId",
        "startNumberOfDaysSinceDiagnosis",
        "endNumberOfDaysSinceDiagnosis",
        "eventType",
    ]
    current_timeline_attr: dict = {
        "Clinical Status": ["CLINICAL_EVENT_TYPE"],
        "IMAGING": [
            "DIAGNOSTIC_TYPE",
            "DIAGNOSTIC_TYPE_DETAILED",
            "BODY_PART",
            "FLYWHEEL_URL",
        ],
        "SPECIMEN": ["SAMPLE_ID", "SPECIMEN_SITE", "SPECIMEN_TYPE"],
        "SURGERY": ["SUBTYPE", "EXTENT_OF_TUMOR_RESECTION"],
        "TREATMENT": ["TREATMENT_TYPE", "SUBTYPE", "AGENT"],
    }

    current_timeline_data: dict = {}
    for attr in current_timeline_attr:
        current_timeline_data[attr] = []
    study_timeline_data: list = requests.get(
        f"{url}/api/studies/{study_id}/clinical-events?projection=DETAILED",
        headers=auth_headers,
        timeout=360,
    ).json()
    for entry in study_timeline_data:
        event_type: str = entry["eventType"]
        temp_event: list = []
        for common in common_fields:
            if common in entry:
                temp_event.append(entry[common])
            else:
                temp_event.append("")
        # event-specific attributes have funky organization,like 'attributes': 
        # [{'key': 'CLINICAL_EVENT_TYPE', 'value': 'Initial CNS Tumor'}],
        # easier to flatten it into a normal dict:
        flatten_attrs: dict = {x["key"]: x["value"] for x in entry["attributes"]}
        event_specific: list = current_timeline_attr[event_type]
        for field in event_specific:
            if field in flatten_attrs:
                temp_event.append(flatten_attrs[field])
            else:
                temp_event.append("")
        current_timeline_data[event_type].append("\t".join(map(str, temp_event)))
    return current_timeline_data


def output_file_by_id(
    select_id: set,
    update_list: list,
    header: list[str],
    id_field: str,
    outfile_name: str,
    big_head: list | None = None,
) -> None:
    """Given a list containing tab-separated entries, select those entries where a field matches a key from a set and print those selections to a file.

    select_id: set of IDs to filter down update_list for output
    update_list: list of entries from update flat file
    header: list containing column name strings
    outfile_name: str output file name
    big_head: str lead header lines, if applicable
    """
    id_index: int = header.index(id_field)
    update_list = [item for item in update_list if item.split("\t")[id_index] in select_id]
    with open(outfile_name, "w") as event_delta:
        event_delta.write("\t".join(header) + "\n" if big_head is None else "".join(big_head))
        event_delta.write("".join(update_list))


def compare_timeline_data(
    current_timeline: dict,
    update_timeline: dict,
    out_dir: str,
    header_info: dict,
    file_ext: dict,
) -> None:
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
        # set to empty to match behavior of current_timeline dict generation
        if event not in update_timeline:
            update_timeline[event] = []
        # strip new line
        update_ids = {item.rstrip("\n") for item in update_timeline[event]}
        diff_ids: set = update_ids - current_ids
        if len(diff_ids) > 1:
            print(
                f"{len(diff_ids)} changes in {event} found for this study. Outputting delta files"
            )
            suffix: str = event_ext_dict[event]
            outfilename: str = f"{out_dir}/data_clinical_timeline_{suffix}"
            uniq_patients: set[str] = {item.split("\t")[0] for item in diff_ids}
            output_file_by_id(
                select_id=uniq_patients,
                update_list=update_timeline[event],
                header=header_info[event],
                id_field="PATIENT_ID",
                outfile_name=outfilename,
                big_head=None,
            )


def run_py(args):
    # Create incremental study updates dirs
    delta_dir: str = f"{args.study}_delta_data"
    os.makedirs(delta_dir, exist_ok=True)
    add_dir: str = f"{args.study}_add_data/datasheets"
    os.makedirs(add_dir, exist_ok=True)

    # Set Up URL pulling
    print("Processing data clinical info", file=sys.stderr)
    with open(args.token, "r") as token_file:
        token: str = token_file.read().rstrip().split(": ")[1]
    auth_headers: dict = {"Authorization": f"Bearer {token}"}
    url_object = urlparse(args.url)
    formatted_url: str = f"{url_object.scheme}://{url_object.hostname}"

    # Fetch attribute keys
    fetch_attr_path: str = "api/clinical-attributes/fetch"
    attr_keys: list[dict[str, str]] = requests.post(
        f"{formatted_url}/{fetch_attr_path}",
        headers=auth_headers,
        json=[args.study],
        params={"projection": "ID"},
        timeout=30,
    ).json()

    # Hardcoding a JSON-like payload that will dictate do similar comparisons
    # For each comparison it has implicit and skip attributes
    # Also includes the name of the data table that has the relevant information as well as the key in that table to locate the information
    aggregate_vals: list[str] = ["SPECIMEN_ID", "EXPERIMENT_STRATEGY"]
    # get all datasheets
    datasheet_dir: str = args.datasheets.rstrip("/")
    class Data(TypedDict):
        attr_implicit: list[str]
        attr_skip: list[str]
        data_table: str
        data_table_key: str
    comparisons: dict[str, Data] = {
        "SAMPLE": {
            "attr_implicit": ["PATIENT_ID"],
            "attr_skip": ["FRACTION_GENOME_ALTERED", "MUTATION_COUNT"],
            "data_table": f"{datasheet_dir}/data_clinical_sample.txt",
            "data_table_key": "SAMPLE_ID",
        },
        "PATIENT": {
            "attr_implicit": [],
            "attr_skip": ["SAMPLE_COUNT"],
            "data_table": f"{datasheet_dir}/data_clinical_patient.txt",
            "data_table_key": "PATIENT_ID",
        },
    }

    for comp in comparisons:
        if comp == "PATIENT":
            current_attr_keys: set[str] = {
                x["clinicalAttributeId"] for x in attr_keys if x["patientAttribute"]
            }
        else:
            current_attr_keys: set[str] = {
                x["clinicalAttributeId"] for x in attr_keys if not x["patientAttribute"]
            }

        # Remove skip attrs and add implicit attrs
        current_attr_keys -= set(comparisons[comp]["attr_skip"])
        current_attr_keys |= set(comparisons[comp]["attr_implicit"])

        # Pull current data from portal URL
        current_data: dict = data_clinical_from_study(
            url=formatted_url,
            auth_headers=auth_headers,
            study_id=args.study,
            data_type=comp,
            aggr_list=aggregate_vals,
        )

        # build update data from file
        update_header, update_columns, update_body = parse_file(
            file_path=comparisons[comp]["data_table"], header_symbol="#", column_names=True,
        )

        # Build update_attr_keys from column names
        update_attr_keys: set[str] = set(update_columns)
        update_attr_keys.remove(comparisons[comp]["data_table_key"])

        # Get attrs uniq to the current and update datasets
        current_only_attr: set[str]
        update_only_attr: set[str]
        shared_attr: set[str]
        current_only_attr, update_only_attr, shared_attr = get_set_venn(
            left=current_attr_keys, right=update_attr_keys
        )

        delta_data, delta_lines, update_only_data, update_only_lines = split_data_file(
            data_lines=update_body,
            column_names=update_columns,
            key=comparisons[comp]["data_table_key"],
            aggr_list=aggregate_vals,
            current_data=current_data,
            shared_attrs=shared_attr,
        )

        # Get ids uniq to the current and update datasets
        update_only_ids: set[str] = set(update_only_data.keys())
        shared_ids: set[str] = set(shared_data.keys())
        current_only_ids: set[str] = set(current_data.keys()) - shared_ids

        # report before and after for deltas
        delta_attr_count: dict[str, int] = {}
        if delta_data:
            delta_attr_count = print_and_count_delta_attr(
                clin_type=clin_type,
                delta_data=delta_data,
                current_data=current_data,
                shared_attr=shared_attr,
                out_filename=f"{comp.lower()}_current_v_update.txt",
            )
            print_parsed_file(
                header=update_header,
                body=delta_lines,
                out_filename=f"{delta_dir}/data_clinical_{comp.lower()}.txt"
            )

        if update_only_data:
            print_parsed_file(
                header=update_header,
                body=update_only_lines,
                out_filename=f"{add_dir}/data_clinical_{comp.lower()}.txt"
            )
            if comp == "SAMPLE":
                _, etl_header, etl_data = parse_file(
                    file_path=args.manifest, header_symbol=None, column_names=True
                )
                output_file_by_id(
                    select_id=update_only_ids,
                    update_list=etl_data,
                    header=etl_header,
                    id_field="Cbio_Tumor_Name",
                    outfile_name=f"{args.study}_add_data/cbio_file_name_id.txt",
                )

        output_clinical_changes(
            clin_type=comp,
            current_only_ids=current_only_ids,
            update_only_ids=update_only_ids,
            current_attr_only=current_only_attr,
            update_attr_only=update_only_attr,
            attr_cts=delta_attr_cts,
            suffix=f"{comp}.txt",
        )

    # timeline-diffs
    timeline_update_files: list = glob.glob(f"{datasheet_dir}/data_clinical*timeline*")
    if len(timeline_update_files):
        print("Getting current timeline data, please wait...", file=sys.stderr)
        timeline_event_dict: dict = {
            "clinical_event.txt": "Clinical Status",
            "imaging.txt": "IMAGING",
            "specimen.txt": "SPECIMEN",
            "surgery.txt": "SURGERY",
            "treatment.txt": "TREATMENT",
        }
        timeline_current: dict = data_clinical_timeline_from_current(
            url=formatted_url, auth_headers=auth_headers, study_id=args.study
        )
        timeline_update_file_prefix: str = f"{datasheet_dir}/data_clinical_timeline_"
        # store file data in memory by event type to match current pull dict created
        timeline_update: dict = {}
        # store file headers to repurpose
        event_file_head: dict = {}
        for path in timeline_update_files:
            suffix: str = path.replace(timeline_update_file_prefix, "")
            event_type: str = timeline_event_dict[suffix]
            _, event_file_head[event_type], timeline_update[event_type] = parse_file(
                file_path=path, header_symbol=None, column_names=True
            )
        print("Comparing current timeline to update", file=sys.stderr)
        compare_timeline_data(
            current_timeline=timeline_current,
            update_timeline=timeline_update,
            out_dir=delta_dir,
            header_info=event_file_head,
            file_ext=timeline_event_dict,
        )

    print("Done!", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Compare update clinical and timeline data to current cbioportal instance. \
        Outputs changes summaries, id lists, and delta files. \
        Recommend running cbioportal_etl/scripts/get_study_metadata.py to get file inputs",
    )
    parser.add_argument(
        "-u",
        "--url",
        action="store",
        dest="url",
        help="url to search against",
        default="https://pedcbioportal.kidsfirstdrc.org/api/v2/api-docs",
    )
    parser.add_argument(
        "-s",
        "--study",
        action="store",
        dest="study",
        help="Cancer study ID to compare on server",
    )
    parser.add_argument(
        "-t",
        "--token",
        action="store",
        dest="token",
        help="Token file obtained from Web API",
    )
    parser.add_argument(
        "-ds",
        "--datasheet-dir",
        action="store",
        dest="datasheets",
        help="Directory containing cBio-formatted  metadata, typically named data_clinical_*",
    )
    parser.add_argument(
        "-m",
        "--manifest",
        default="cbio_file_name_id.txt",
        help="Manifest file (default: cbio_file_name_id.txt from Step 1 output)",
    )
    args = parser.parse_args()
    run_py(args)


if __name__ == "__main__":
    main()
