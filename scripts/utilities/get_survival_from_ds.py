#!/usr/bin/env python3

import sys
from requests import request
import argparse
import concurrent.futures
import json
from time import sleep

def query_ds(pt_id):
    try:
        pt_id = pt_id.rstrip('\n')
        outcome_url = url + "outcomes?participant_id=" + pt_id
        out_info = request('GET', outcome_url)
        statuses = []
        for result in out_info.json()['results']:
            statuses.append([pt_id, result['vital_status'], str(result['age_at_event_days'])])
        return statuses
    except Exception as e:
        sys.stderr.write(str(e))
        sys.stderr.write('Could not get outcome info for ' + pt_id + "\n")
        sys.exit(1)


parser = argparse.ArgumentParser(description='Script to walk through data service and grab all relevant meetadata'
                                             ' by bs id.')
parser.add_argument('-u', '--kf-url', action='store', dest='url',
                    help='Kids First data service url, i.e. https://kf-api-dataservice.kidsfirstdrc.org/')
parser.add_argument('-p', '--pt-list', action='store', dest='pt_list',
                    help='pt_id list')

args = parser.parse_args()

url = args.url
pt_file = open(args.pt_list)
out = open('outcomes.txt', 'w')

x = 1
m = 100

with concurrent.futures.ThreadPoolExecutor(16) as executor:
    results = {executor.submit(query_ds, entry): entry for entry in pt_file}
    for result in concurrent.futures.as_completed(results):
        if x % m == 0:
            sys.stderr.write('Processed ' + str(x) + ' lines\n')
            sys.stderr.flush()
        outcomes = result.result()
        if len(outcomes) > 0:
            for outcome in outcomes:
                out.write("\t".join(outcome) + "\n")
        else:
            sys.stderr.write('No outcome for a patient\n')
        x += 1
