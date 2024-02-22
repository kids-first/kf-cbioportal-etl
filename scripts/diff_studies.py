"""
Script to check a study on pedcbioportal for differences against a local build
"""

import argparse
from bravado.client import SwaggerClient
from bravado.requests_client import RequestsClient
from urllib.parse import urlparse
import pdb


def main():
    parser = argparse.ArgumentParser(
        description="Compare local clinical data to server"
    )
    parser.add_argument(
        "-u",
        "--url",
        action="store",
        dest="url",
        help="url to search against",
        default="https://pedcbioportal.kidsfirstdrc.org/api/v2/api-docs"
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
    # object is a big ass array, one entry per attribute, per patient/sample? - so like if a table were just concatenated into a single vector 
    clinical_data = cbioportal.Clinical_Data.getAllClinicalDataInStudyUsingGET(studyId=study).result()
    # would need to loop through patient IDs get patient-specific fields, like age?
    # get_clin_pt = cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(studyId=study,patientId='PT_00G007DM' ).result()
    pdb.set_trace()
    hold=1

if __name__ == '__main__':
    main()
