#!/usr/bin/env python3

import sevenbridges as sbg
import sys
import argparse
import concurrent.futures
from sevenbridges.http.error_handlers import rate_limit_sleeper, maintenance_sleeper
import re
import json


def process_project_tasks(project, task, atype, cna_flag=None):
    task_name = task.name.rstrip('\n')
    task_id = task.id
    if atype == 'DNA' and 'strelka2_vep_maf' in task.outputs and task.status == 'COMPLETED' and \
            re.search(config_data['dna_task_keyword'], task.name):
        try:
            (bs_id1, bs_id2) = re.findall('BS_[A-Z0-9]*',task.name)
            tum_bs_id = bs_id1
            norm_bs_id = bs_id2
            if config_data['tn_order'] == 'NT':
                tum_bs_id = bs_id2
                norm_bs_id = bs_id1
            outs = task.outputs['strelka2_vep_maf'].name
            if cna_flag == 1:
                outs += ',' + task.outputs['cnv_pval'].name
            return '\t'.join((tum_bs_id, norm_bs_id, task_name, task_id, atype, outs, project)) + '\n'
        except Exception as e:
            sys.stderr.write(str(e) + '\n')
            sys.stderr.write('Failed to process ' + task.id + ' ' + task.name + ', skipping\n')
    elif atype == 'RNA':
        try:
            if task.outputs is not None and 'RSEM_gene' in task.outputs:
                bs_id = task.outputs['RSEM_gene'].metadata['Kids First Biospecimen ID']
                outs = task.outputs['RSEM_gene'].name
                return '\t'.join((bs_id, 'NA', task_name, task_id, atype, outs, project)) + '\n'

        except Exception as e:
            sys.stderr.write(str(e) + '\n')
            sys.stderr.write('Task ' + task_name + ' with ID ' + task_id + ' failed to gather all info\n')
            sys.stderr.flush()
            return None


def get_file_info(api, project, atype, out_fh, cna_flag=None):
    proj_tasks = api.tasks.query(project=project).all()
    sys.stderr.write('Processing project ' + project + '\n')
    x = 1
    n = 20
    with concurrent.futures.ThreadPoolExecutor(config_data['cpus']) as executor:
        results = {executor.submit(process_project_tasks, project, task, atype, cna_flag): task for task in proj_tasks}
        for result in concurrent.futures.as_completed(results):
            if x % n == 0:
                sys.stderr.write('Processed ' + str(x) + ' tasks\n')
                sys.stderr.flush()
            outstr = result.result()
            if outstr is not None:
                out_fh.write(outstr)
            x += 1


parser = argparse.ArgumentParser(description='Get all relevant analyzed file outputs from projects in cavatica.')
parser.add_argument('-o', '--output', action='store', dest='out', help='output file name')
parser.add_argument('-p', '--profile', action='store', dest='profile', help='cavatica profile name. requires '
                                                                            '.sevenbridges/credentials file be present')
parser.add_argument('-c', '--config', action='store', dest='config_file', help='json config file with data types and '
                                                                               'data locations')

args = parser.parse_args()
config = sbg.Config(profile=args.profile)
api = sbg.Api(config=config, error_handlers=[
                sbg.http.error_handlers.rate_limit_sleeper,
                sbg.http.error_handlers.maintenance_sleeper])
with open(args.config_file) as f:
    config_data = json.load(f)
dna_types = config_data['cav_dna_projects']
rna_types = config_data['cav_rna_projects']
out_fh = open(args.out, 'w')
out_fh.write('T/CL BS ID\tNorm BS ID\tTask Name\tTask ID\tAnalyte Type\tRelevant Outputs\tSource Project\n')
for project in dna_types:
    get_file_info(api, project, 'DNA', out_fh, config_data['cna_flag'])
for project in rna_types:
    get_file_info(api, project, 'RNA', out_fh)
out_fh.close()
sys.stderr.write('Complete\n')

