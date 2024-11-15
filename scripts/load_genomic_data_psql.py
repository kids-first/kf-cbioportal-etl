#!/usr/bin/env python3
"""
This script takes a cBio-formatted genomic file and loads it into a database table.
You can either initialize a table using "create" mode, add data/samples using "append".
YOU MUST CREATE YOUR OWN .ini FILE WITH LOGIN PARAMS BEFORE RUNNING
"""
import psycopg
from psycopg import sql
from get_study_metadata import config
import argparse
import pdb
import sys


def format_tbl_name(tbl_path):
    """
    Formats the table name to fit SQL standard of psycopg package
    """

    schema = None
    tname = '"tbl_path"'
    formatted_tpath = sql.Identifier(tbl_path)
    if "." in tbl_path:
        (schema, tname) = tbl_path.split(".")
        formatted_tpath =sql.Identifier(schema, tname)
        tname = f"""{schema}."{tname}" """
    return formatted_tpath, tname


def mt_execute_data_load(line, formatted_tpath, header, db_conn):
    try:
        info = line.rstrip('\n').split('\t')
        vals = [float(old) for old in info[1:]]
        db_conn.execute(sql.SQL("INSERT INTO {name} ({primary_key}, samples, data) VALUES (%s, %s, %s)").format(name=formatted_tpath, primary_key=sql.Identifier(header[0])), [info[0], header[1:], vals])
        return 0
    except psycopg.Error as e:
        print(e, file=sys.stderr)
        print("Failed data load for primary key {}".format(info[0]), file=sys.stderr)
        sys.stderr.flush()
        db_conn.close()
        return 1


def create_table(header, data, tbl_path, data_type, formats, db_conn, wtypes):
    """
    This function creates the desired table from a formatted cBio file depending on it's table format:
    wide (gene by sample name) or long (tsv with headers) then loads accordingly 
    """
    if formats[data_type] == "wide":
        print("Creating {} type table for {} type data, starting with table def".format(formats[data_type], data_type), file=sys.stderr)
        (formatted_tpath, tname) = format_tbl_name(tbl_path)
        print("DEBUG: Waiting for pool connection", file=sys.stderr)
        create_sql = sql.SQL("CREATE TABLE IF NOT EXISTS {name} ({primary_key} PRIMARY KEY, samples varchar[], {data})").format(name=formatted_tpath, primary_key=sql.SQL(" {} {}").format(sql.Identifier(header[0]), sql.SQL("varchar")),data=sql.SQL(" {} {}").format(sql.Identifier("data"), sql.SQL(wtypes[data_type])))
        db_conn.execute(create_sql)
        num_records = 1
        batch_size = 2000
        print("Loading data from TSV in batches of {} records".format(batch_size), file=sys.stderr)
        batch_load = []
        load_sql = f"""INSERT INTO {tname} VALUES (%s, %s, %s);"""
        for line in data:
            if num_records % batch_size == 0:
                print("Read {} records into memory".format(batch_size))
                db_conn.cursor().executemany(load_sql, batch_load)
                print("Loaded {} records into table".format(num_records))
                batch_load = []
            info = line.rstrip('\n').split('\t')
            batch_load.append( (info[0], header[1:], [float(old) for old in info[1:]] ) )
            num_records += 1
        if num_records % batch_size != 0:
            db_conn.cursor().executemany(load_sql, batch_load)
        print("Batch load into table completed!", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Pull clinical data and genomics file etl support from D3b data warehouse.")

    parser.add_argument("-d", "--db-ini", action="store", dest="db_ini", help="Database config file - formatting like aws or sbg creds")
    parser.add_argument("-p", "--profile", action="store", dest="profile", help="ini profile name", default="postgresql")
    parser.add_argument("-m", "--mode", action="store", dest="mode", help="Load mode - create or append", default="create")
    parser.add_argument("-n", "--table-name", action="store", dest="name", help="Name of table, including schema")
    parser.add_argument("-i", "--input", action="store", dest="input", help="Input tsv file")
    parser.add_argument("-t", "--type", action="store", dest="type", help="cBio file type, i.e. cnv, maf, rsem, sv")

    args = parser.parse_args()

    formats = { "cnv": "wide", "rsem": "wide", "maf": "long", "sv": "long" }
    wtypes = { "cnv": "integer[]", "rsem": "double precision[]" }

    with open(args.input) as file_in:
        head = next(file_in)
        header = head.rstrip('\n').split('\t')

        # Load database login info
        db = config(filename=args.db_ini, section=args.profile)
        psql_uri = "postgresql://{user}:{password}@{host}/{dbname}".format(user=db['user'], password=db['password'], host=db['host'], dbname=db['database'])
        print("Establishing connection", file=sys.stderr)
        with psycopg.connect(conninfo=psql_uri, autocommit=True) as db_conn:
            if args.mode == 'create':
                create_table(header, file_in, args.name, args.type, formats, db_conn, wtypes)
                print("Table {} created!".format(args.name))
            elif args.mode == 'append':
                # append_table(header, file_in, args.name, args.type, formats, cur, conn, wtypes)
                print("Append to {} completed!".format(args.name))
            else:
                print("ERROR, mode choices are create or append, you chose {mode}".format(mode=args.mode))
                exit(1)

if __name__ == "__main__":
    main()
