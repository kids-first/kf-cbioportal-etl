"""
Script to check a study on pedcbioportal for differences against a local build
"""

import argparse
from bravado.client import SwaggerClient
from bravado.requests_client import RequestsClient
from urllib.parse import urlparse
import pdb


def clinical_diffs(portal, build, portal_attr, build_attr, clin_type, out):
    """
    Compare differences in portal sample data and build.
    """
    # gross  ID diffs
    portal_clinical_ids = set(list(portal.keys()).sort())
    build_clinical_ids = set(list(build.keys()).sort())
    portal_only = list(portal_clinical_ids - build_clinical_ids)
    build_only = list(build_clinical_ids - portal_clinical_ids)
    common_samp_ids = list(portal_clinical_ids & build_clinical_ids)
    # gross attribute diffs
    portal_attr_only = list(portal_attr - build_attr)
    build_attr_only = list(build_attr - portal_attr)
    common_attr = list(portal_attr & build_attr)

    # focus on common samp and common attr, as "everything is different for x" is not that useful
    print("Per " + clin_type + " changes:", file=out)
    attr_cts = {}
    for samp_id in common_samp_ids:
        for attr in common_attr:
            if portal[samp_id][attr] != build[samp_id][attr]:
                print("{} {} attribute {} would change from {} to {}".format(clin_type, samp_id, attr, portal[samp_id][attr], build[samp_id][attr]), file=out)
                if attr not in attr_cts:
                    attr_cts[attr] = 0
                attr_cts[attr] += 1
    print("CHANGE SUMMARY:", file=out)
    if len(portal_only) > 0:
        print("{} {}s in portal would be removed: {}".format(len(portal_only), clin_type, ",".join(portal_only)), file=out)
    if len(build_only) > 0:
        print("{} {}s in build would be added to the portal: {}".format(len(build_only), clin_type,  ",".join(build_only)), file=out)
    if len(portal_attr_only) > 0:
        print("{} attributes in portal would be removed: {}".format(len(portal_attr_only), ",".join(portal_attr_only)), file=out)
    if len(build_attr_only) > 0:
        print("{} attributes in build would be added to the portal: {}".format(len(build_attr_only), ",".join(build_attr_only)), file=out)
    for attr in attr_cts:
        print("{} has {} change(s)".format(attr, attr_cts[attr]), file=out)


def split_sort_field(value, sep):
    """
    For ease of comparison, aggregate attributes concatenated by a separator may not be in the same order, but order does not matter.
    Therefore, sort them so that when compared, no errors are triggered
    """
    return sep.join(value.split(sep).sort())


def table_to_dict(in_file, key, aggr_list):
    """
    Take a text file and convert to dict with certain row value as primary, all other row values as subkeys.
    Also return a set of attribute keys
    """
    with open(in_file) as f:
        # skip lines starting with hash until normal header is reached
        for entry in f:
            if entry[0] != "#":
                head = next(f)
                header = head.rstrip('\n').split('\t')
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
            data_dict[data[primary]] = {}
            # sort aggr fields
            for i in aggr_head:
                entry[i] = split_sort_field(entry[i], ";")
            # two loops, for up until primary key, then after
            for i in range(0, primary, 1):
                data_dict[data[primary]][header[i]] = entry[i]
            for i in range((primary + 1), len(entry), 1):
                data_dict[data[primary]][header[i]] = entry[i]
    return data_dict, set(header.sort())


def data_clinical_from_study(cbio_conn, study_id, data_type, aggr_list):
    """
    Get all the column-value pairs for each data_type(SAMPLE or PATIENT) for a specific study
    Convert result to dict
    Also return a set of attribute keys
    """
    # object is a big ass array of struct, one entry per attribute, per patient/sample - so like if a table were just concatenated into a single vector 
    data_clinical = cbio_conn.Clinical_Data.getAllClinicalDataInStudyUsingGET(studyId=study_id, projection='DETAILED', clinicalDataType=data_type).result()
    data_dict = {}
    attr_keys = []
    # Use sampleId or patientID, clinicalAttributeId (column name) and value
    attr_dict = {"SAMPLE": "sampleId", "PATIENT": "patientId" }
    for entry in data_clinical:
        clinical_id = getattr(entry, attr_dict[data_type])
        if clinical_id not in data_dict:
            data_dict[clinical_id] = {}
        value = entry.value
        attr_id = entry.clinicalAttributeId
        if attr_id in aggr_list:
            value = split_sort_field(value, ";")
        data_dict[clinical_id][attr_id] = value
        # need to fix this...can probably just use api to get attribute keys instead
        attr_keys.append(attr_id)
    return data_dict, set(attr_keys.sort())


def main():
    parser = argparse.ArgumentParser(
        description="Compare local clinical data to server"
    )
    parser.add_argument(
        "-u", "--url", action="store", dest="url", help="url to search against", default="https://pedcbioportal.kidsfirstdrc.org/api/v2/api-docs"
    )
    parser.add_argument(
        "-s", "--study", action="store", dest="study", help="Cancer study ID to compare on server"
    )
    parser.add_argument(
        "-t", "--token", action="store", dest="token", help="Token file obtained from Web API"
    )
    parser.add_argument(
        "-d", "--datasheet-dir", action="store", dest="data_dir", help="Directory containing data_clinical_*.txt"
    )

    args = parser.parse_args()
    with open(args.token, 'r') as token_file:
        token = token_file.read().rstrip().split(': ')[1]

    url_object = urlparse(args.url)

    http_client = RequestsClient()
    http_client.set_api_key(
        '{}'.format(url_object.hostname), 'Bearer {}'.format(token),
        param_name='Authorization', param_in='header'
    )

    cbioportal = SwaggerClient.from_url(args.url,
                                        http_client=http_client,
                                        config={"validate_requests":False,
                                                "validate_responses":False,
                                                "validate_swagger_spec": False}
    )
    
    portal_sample_data, portal_sample_attr_keys = data_clinical_from_study(cbioportal, args.study, "SAMPLE")
    build_sample_data, build_sample_attr_keys = table_to_dict(args.data_dir + "/data_clinical_sample.txt")
    sample_diff_out = open('sample_portal_v_build.txt', 'w')
    clinical_diffs(portal_sample_data, build_sample_data, portal_sample_attr_keys, build_sample_attr_keys, "Sample", sample_diff_out)
    sample_diff_out.close()
    # hardcode for now names of aggregate fields
    aggr_list = ["SPECIMEN_ID", "EXPERIMENT_STRATEGY"]
    portal_patient_data, portal_patient_attr_keys =  data_clinical_from_study(cbioportal, args.study, "PATIENT", aggr_list)
    build_patient_data, build_patient_attr_keys = table_to_dict(args.data_dir + "/data_clinical_patient.txt", aggr_list)
    patient_diff_out = open('patient_portal_v_build.txt', 'w')
    clinical_diffs(portal_patient_data, build_patient_data, portal_patient_attr_keys, build_patient_attr_keys, "Patient", patient_diff_out)
    patient_diff_out.close()

    pdb.set_trace()
    hold=1

if __name__ == '__main__':
    main()
