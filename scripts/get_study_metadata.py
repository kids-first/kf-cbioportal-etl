"""
This script downloads all relevant data_clinical and genomics file etl support files from the D3b Data Warehouse.
YOU MUST CREATE YOUR OWN .ini FILE WITH LOGIN PARAMS BEFORE RUNNING
"""
import psycopg2
from psycopg2 import sql
import os
from configparser import ConfigParser
import argparse
import json
import pdb


def config(filename='database.ini', section='postgresql'):
    # stolen from https://www.postgresqltutorial.com/postgresql-python/connect/
    # create a parser
    parser = ConfigParser()
    # read config file
    parser.read(filename)

    # get section, default to postgresql
    db = {}
    if parser.has_section(section):
        params = parser.items(section)
        for param in params:
            db[param[0]] = param[1]
    else:
        raise Exception('Section {0} not found in the {1} file'.format(section, filename))

    return db


def generic_pull(db_cur, tbl_name):
    """
    Simple helper function for SELECT * database calls
    """
    if '.' not in tbl_name:
        tbl_sql = sql.SQL('SELECT * FROM {};').format(sql.Identifier(tbl_name))
    else:
        (schema, table) = tbl_name.split('.')
        tbl_sql = sql.SQL('SELECT * FROM {}.{};').format(sql.Identifier(schema), sql.Identifier(table))

    db_cur.execute(tbl_sql)
    rows = db_cur.fetchall()
    colnames = [desc[0] for desc in db_cur.description]
    return (rows, colnames)


def generic_print(out_file, rows, colnames):
    """
    Simple helper function to print results to existing file handle as tsv
    """
    out_file.write("\t".join(colnames) + "\n")
    for row in rows:
        # convert None to empty str
        new_row = [str(i or '') for i in row]
        out_file.write("\t".join(new_row) + "\n")
    out_file.close()
    return 0


def get_data_clinical(db_cur, config_dict, prefix):
    """
    Depending on the prefix of patient or sample, will pull from related tables,
    only use related header info present in table, and print the combined results.
    This is because cBioportal flat file had 5 header rows and it can vary between projects which are used
    """
    if prefix not in ['sample', 'patient']:
        print("Please pass sample or patient for prefix.")
        exit(1)
    tbl_name = config_dict['database_pulls'][prefix + '_file']['table']

    # get table contents
    (rows, colnames) = generic_pull(db_cur, tbl_name)

    # use table header from colnames, and use to select file header
    head_name = config_dict['database_pulls'][prefix + '_head']['table']
    # get sample table contents, have to split if format schema.table
    if '.' not in head_name:
        head_sql = sql.SQL('SELECT {} FROM {};').format(sql.SQL(',').join(map(sql.Identifier, colnames)), sql.Identifier(head_name))
    else:
        (schema, table) = head_name.split('.')
        head_sql = sql.SQL('SELECT {} FROM {}.{};').format(sql.SQL(',').join(map(sql.Identifier, colnames)), sql.Identifier(schema), sql.Identifier(table))
    db_cur.execute(head_sql)
    head = db_cur.fetchall()
    # create output file and combine results for final product
    out_file = open(datasheet_dir + "/" + config_data['database_pulls'][prefix + '_file']['out_file'], 'w')
    for row in head:
        out_file.write("\t".join(row) + "\n")
    generic_print(out_file, rows, colnames)
    return 0


def get_manifests(db_cur, config_dict):
    """
    This iterates through the manifests section of the database_pulls and grabs and outputs all listed file manifests for ec2 download
    """
    # Just a pointer variable for ease of reading and coding
    manifests = config_dict['database_pulls']['manifests']
    for manifest in manifests:
        try:
            tbl_name = manifests[manifest]['table']
            file_types = manifests[manifest]['file_type']
            if '.' not in tbl_name:
                manifest_sql = sql.SQL('SELECT * FROM {} WHERE file_type in ({});').format(sql.Identifier(tbl_name), sql.SQL(',').join(map(sql.Literal, file_types)))
            else:
                (schema, table) = tbl_name.split('.')
                manifest_sql = sql.SQL('SELECT * FROM {}.{} WHERE file_type in ({});').format(sql.Identifier(schema), sql.Identifier(table), sql.SQL(',').join(map(sql.Literal, file_types)))
            db_cur.execute(manifest_sql)
            rows = db_cur.fetchall()
            colnames = [desc[0] for desc in db_cur.description]
            out_file = open(manifests[manifest]['out_file'], 'w')
            generic_print(out_file, rows, colnames)

        except Exception as e:
            print(e)
            exit(1)
    return 0


parser = argparse.ArgumentParser(description="Pull clinical data and genomics file etl support from D3b data warehouse.")

parser.add_argument("-d", "--db-ini", action="store", dest="db_ini", help="Database config file - formatting like aws or sbg creds")
parser.add_argument("-p", "--profile", action="store", dest="profile", help="ini profile name", default="postgresql")
parser.add_argument("-c", "--config", action="store", dest="config_file", help="json config file with meta information; see REFS/pbta_all_case_meta_config.json example",)

args = parser.parse_args()
# Load database login info
params = config(filename=args.db_ini, section=args.profile)
datasheet_dir = 'datasheets'
# Load json config file with database pull info
with open(args.config_file) as f:
    config_data = json.load(f)
try:
    conn = psycopg2.connect(**params)
    cur = conn.cursor()
        
    # dict to track keys with specific database calls
    special_keys = {"sample_head": 0, "sample_file": 0, "patient_head": 0, "patient_file": 0, "manifests": 0}

    try:
        os.mkdir(datasheet_dir)
    except Exception as e:
        print(str(e) + ' IGNORE!')
    get_data_clinical(cur, config_data, 'sample')
    get_data_clinical(cur, config_data, 'patient')
    get_manifests(cur, config_data)

    # For all other tables to be printed simply, not in special_keys
    for key in config_data['database_pulls']:
        if key not in special_keys and key != '_comment':
            tbl_name = config_data['database_pulls'][key]['table']
            (rows, colnames) = generic_pull(cur, tbl_name)
            out_fn = config_data['database_pulls'][key]['out_file']
            if key == 'gene_file':
                out_fn = datasheet_dir + "/" + out_fn
            out_file = open(out_fn, 'w')
            generic_print(out_file, rows, colnames)

    # close the communication with the PostgreSQL
    cur.close()
except (Exception, psycopg2.DatabaseError) as error:
    print(error)
    conn.close()
conn.close()