#!/usr/bin/env python3
"""
This script takes a cBio-formatted genomic file and loads it into a database table.
You can either initialize a table using "create" mode, add data/samples using "append".
YOU MUST CREATE YOUR OWN .ini FILE WITH LOGIN PARAMS BEFORE RUNNING
"""
import psycopg2
from psycopg2 import sql
import argparse
import pdb
import sys
from get_study_metadata import config


def process_wide_header(header, tbl_path, wtypes, data_type):
    """
    This function takes the input file header, path to table, data type dict and data type.
    In return, for downstream use, it gives the schema (if any), table name, sql-object-formatted table path and array of column-sql data type tuples
    """
    col_type = []
    col_type.append((header[0], "varchar PRIMARY KEY"))
    for cname in header[1:]:
        col_type.append((cname, wtypes[data_type]))
    schema = None
    tname = tbl_path
    formatted_tpath = sql.Identifier(tbl_path)
    if "." in tbl_path:
        (schema, tname) = tbl_path.split(".")
        formatted_tpath =sql.Identifier(schema, tname)
    return schema, tname, formatted_tpath, col_type


def create_table(header, data, tbl_path, data_type, formats, db_cur, db_conn, wtypes):
    """
    This function creates the desired table from a formatted cBio file depending on it's table format:
    wide (gene by sample name) or long (tsv with headers) then loads accordingly 
    """
    if formats[data_type] == "wide":
        # Used DuckAssist (AI) to help write this
        print("Creating {} type table for {} type data, starting with table def".format(formats[data_type], data_type), file=sys.stderr)
        (schema, tname, formatted_tpath, col_type) = process_wide_header(header, tbl_path, wtypes, data_type)
        create_sql = sql.SQL("CREATE TABLE IF NOT EXISTS {name} ({cols})").format(name=formatted_tpath, cols=sql.SQL(", ").join(sql.SQL(" {} {}").format(sql.Identifier(name), sql.SQL(type)) for name, type in col_type ) )
        db_cur.execute(create_sql)
        db_conn.commit()
        print("Loading data from TSV", file=sys.stderr)
        if schema is not None:
            db_cur.execute(f'SET search_path TO {schema}')
        db_cur.copy_from(data, tname, sep='\t')
        db_conn.commit()
        return formatted_tpath, col_type


def append_table(header, data, tbl_path, data_type, formats, db_cur, db_conn, wtypes):
    """
    This function appends the desired data to an existing table table from a formatted cBio file depending on it's table format:
    wide (gene by sample name) or long (tsv with headers) then loads accordingly.
    It creates a temp table using the existing create_table function, add the new sample columns, inserts the new data using the primary key, then removes temp table
    """

    if formats[data_type] == "wide":
        print("Adding data to existing {} format table from {} data".format(formats[data_type], data_type), file=sys.stderr)
        temp_table = "bix_workflows.delme_join_temp"
        print("Creating temp table {}".format(temp_table), file=sys.stderr)
        formatted_tmp_tbl_path, tmp_col_type = create_table(header, data, temp_table, data_type, formats, db_cur, db_conn, wtypes)
        formatted_to_append_path = sql.Identifier(tbl_path)
        if "." in tbl_path:
            (schema, tname) = tbl_path.split(".")
            formatted_to_append_path =sql.Identifier(schema, tname)
        print("Adding new columns to table to be updated", file=sys.stderr)
        add_cols_sql = sql.SQL("ALTER TABLE {insert_tbl} {cols}").format(insert_tbl=formatted_to_append_path, cols=sql.SQL(", ").join(sql.SQL(" ADD COLUMN {} {}").format(sql.Identifier(name), sql.SQL(type)) for name, type in tmp_col_type[1:] ) )
        db_cur.execute(add_cols_sql)
        db_conn.commit()
        print("Updating table with new data", file=sys.stderr)
        update_sql = sql.SQL("UPDATE {insert_tbl} b SET {new_data} FROM {tmp_tbl} a WHERE a.{join_col} = b.{join_col}").format(insert_tbl=formatted_to_append_path, new_data=sql.SQL(", ").join(sql.SQL(" {col_name}=a.{col_name}").format(col_name=sql.Identifier(name_type[0])) for name_type in tmp_col_type[1:] ), tmp_tbl=formatted_tmp_tbl_path, join_col=sql.Identifier("Hugo_Symbol"))
        db_cur.execute(update_sql)
        db_conn.commit()
        print("Removing temp table", file=sys.stderr)
        drop_tmp_sql = sql.SQL("DROP TABLE {tmp_tbl}".format(tmp_tbl=temp_table))
        db_cur.execute(drop_tmp_sql)
        db_conn.commit()
        print("Appending data complete!", file=sys.stderr)


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
            create_table(header, file_in, args.name, args.type, formats, cur, conn, wtypes)
            print("Table {} created!".format(args.name))
        elif args.mode == 'append':
            append_table(header, file_in, args.name, args.type, formats, cur, conn, wtypes)
            print("Append to {} completed!".format(args.name))
        else:
            print("ERROR, mode choices are create or append, you chose {mode}".format(mode=args.mode))
            exit(1)

    conn.close()

if __name__ == "__main__":
    main()
