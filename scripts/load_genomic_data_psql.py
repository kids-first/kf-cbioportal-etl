#!/usr/bin/env python3
"""
This script takes a cBio-formatted genomic file and loads it into a database table
YOU MUST CREATE YOUR OWN .ini FILE WITH LOGIN PARAMS BEFORE RUNNING
"""
import psycopg2
from psycopg2 import sql
import argparse
import sys
from get_study_metadata import config


def create_table(header, data, name, type, formats, db_cur, db_conn):
    if formats[type] == "wide":
        col_type = []
        col_type.append((header[0], "varchar PRIMARY KEY"))
        for cname in header[1:]:
            col_type.append((cname, "integer"))
        tname = name
        tbl_name = sql.Identifier(name)
        if "." in name:
            (schema, tname) = name.split(".")
            tbl_name =sql.Identifier(schema, tname)
        # Used DuckAssist (AI) to help write this
        create_sql = sql.SQL("CREATE TABLE IF NOT EXISTS {name} ({cols})").format(name=tbl_name, cols=sql.SQL(", ").join(sql.SQL(" {} {}").format(sql.Identifier(name), sql.SQL(type)) for name, type in col_type ) )
        db_cur.execute(create_sql)
        db_conn.commit()
        db_cur.execute(f'SET search_path TO {schema}')
        db_cur.copy_from(data, tname, sep='\t')
        db_conn.commit()


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
    wtypes = { "cnv": "integer", "rsem": "double precision" }

    with open(args.input) as file_in:
        head = next(file_in)
        header = head.rstrip('\n').split('\t')

    # Load database login info
        params = config(filename=args.db_ini, section=args.profile) 
        conn = psycopg2.connect(**params)
        cur = conn.cursor()

        if args.mode == 'create':
            create_table(header, file_in, args.name, args.type, formats, cur, conn)

if __name__ == "__main__":
    main()