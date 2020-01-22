#!/usr/bin/env python3
import sys
import os
import argparse
import json
import subprocess

parser = argparse.ArgumentParser(description='Mount sbfs, create soft links bases on file type')
parser.add_argument('-m', '--manifest', action='store', dest='manifest', help='cavatica csv manifest file')
parser.add_argument('-c', '--config', action='store', dest='config_file', help='json config file with data types and '
                                                                               'data locations')
parser.add_argument('-p', '--profile', action='store', dest='profile', help='cavatica profile name. requires '
                                                                            '.sevenbridges/credentials file be present')
args = parser.parse_args()
with open(args.config_file) as f:
    config_data = json.load(f)
os.makedirs('mafs', exist_ok=True)
if config_data['cna_flag'] == 1:
    os.makedirs('cnvs', exist_ok=True)
if config_data['rna_flag'] == 1:
    os.makedirs('rsem', exist_ok=True)

# create mounts and make links
m_dict = {}
manifest = open(args.manifest)
# just last two cols required
head = next(manifest)
cwd = os.getcwd() + "/"
for line in manifest:
    info = line.rstrip('\n').split('\t')
    fnames = info[-2].split(',')
    projects = info[-1]
    project = projects.split(',')
    atype = info[4]
    for i in range(len(project)):
        (u,p) = project[i].split('/')
        if p not in m_dict:
            os.mkdir(p)
            cmd = 'sbfs mount ' + p + ' --project ' + project[i] + ' --profile turbo --read-only'
            subprocess.call(cmd, shell=True)
            mdir = cwd + p + "/projects/" + project[i] + "/"
            m_dict[p] = mdir
        cmd = 'ln -s ' + m_dict[p] + fnames[i]
        if atype == 'DNA':
            if fnames[i]:[-3:] == 'maf':
                cmd += " mafs"
            else:
                cmd += " cnvs"
        else:
            cmd += " rsem"
        subprocess.call(cmd, shell=True)
out = open('sbfs_mount.txt', 'w')
for p in m_dict:
    out.write(m_dict[p] + '\n')
out.close()