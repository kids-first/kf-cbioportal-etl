#!/usr/bin/env python3
"""Downloads all relevant data_clinical and genomics file etl support files from the D3b Data Warehouse.

YOU MUST CREATE YOUR OWN .ini FILE WITH LOGIN PARAMS BEFORE RUNNING
"""

import argparse
import json
import os
import sys
from configparser import ConfigParser
from typing import IO, Any

import psycopg2
from psycopg2 import sql


def config(filename: str = "database.ini", section: str = "postgresql") -> dict[str, str]:
    # stolen from https://www.postgresqltutorial.com/postgresql-python/connect/
    # create a parser
    parser: ConfigParser = ConfigParser()
    # read config file
    parser.read(filename)
    # get section, default to postgresql
    db: dict[str, str] = {}
    if parser.has_section(section):
        params = parser.items(section)
        for param in params:
            db[param[0]] = param[1]
    else:
        err_msg: str = f"Section {section} not found in the {filename} file"
        raise Exception(err_msg)
    return db


def generic_pull(
    db_cur: psycopg2.extensions.cursor, tbl_name: str
) -> tuple[list[tuple[Any]], list[str]]:
    """Format SELECT * database calls."""
    if "." not in tbl_name:
        tbl_sql = sql.SQL("SELECT * FROM {};").format(sql.Identifier(tbl_name))
    else:
        (schema, table) = tbl_name.split(".")
        tbl_sql = sql.SQL("SELECT * FROM {}.{};").format(
            sql.Identifier(schema), sql.Identifier(table)
        )

    db_cur.execute(tbl_sql)
    rows: list[tuple[Any]] = db_cur.fetchall()
    colnames: list[str] = [desc[0] for desc in db_cur.description]
    return (rows, colnames)


def generic_print(out_file: IO, rows: list[tuple[Any]], colnames: list[str]) -> int:
    """Print results to existing file handle as tsv."""
    print("\t".join(colnames), file=out_file)
    for row in rows:
        # convert None to empty str
        new_row: list[str] = ["" if i is None else str(i) for i in row]
        print("\t".join(new_row), file=out_file)
    return 0


def get_data_clinical(
    db_cur: psycopg2.extensions.cursor,
    config_dict: dict,
    prefix: str,
    ref_dir: str,
    datasheet_dir: str,
) -> int:
    """Depending on the prefix of patient or sample, will pull from related tables.

    Only use related header info present in table, and print the combined results.
    This is because cBioportal flat file had 5 header rows and it can vary between projects which are used
    """
    if prefix not in ["sample", "patient"]:
        print("Please pass sample or patient for prefix.", sys.stderr)
        sys.exit(1)
    tbl_name: str = config_dict["database_pulls"][f"{prefix}_file"]["table"]

    # get table contents
    (rows, colnames) = generic_pull(db_cur, tbl_name)

    # use table header from colnames, and use to select file header
    head_file: IO = open(ref_dir + config_dict["database_pulls"][prefix + "_head"]["table"])
    # get and read head file
    head_lines: list[str] = head_file.readlines()
    # create output file and combine results for final product
    out_file = open(
        f"{datasheet_dir}/{config_dict['database_pulls'][f'{prefix}_file']['out_file']}", "w"
    )
    # get indices of matching head lines, then print corresponding cBio header values
    # the last row, and the header of the data clinical table should have overlapping values
    head_search: list[str] = head_lines[-1].rstrip("\n").split("\t")
    col_i: list[int] = [head_search.index(col) for col in colnames]
    for i in range(0, len(head_lines) - 1, 1):
        head = [head_lines[i].rstrip("\n").split("\t")[j] for j in col_i]
        print("\t".join(head), file=out_file)
    generic_print(out_file, rows, colnames)
    out_file.close()
    return 0


def run_py(args):
    # Load database login info
    params: dict[str, str] = config(filename=args.db_ini, section=args.profile)
    study: str = args.study_config
    datasheet_dir: str = "datasheets"
    # Load json config file with database pull info
    with open(study) as f:
        config_data: dict = json.load(f)
    try:
        conn: psycopg2.extensions.connection = psycopg2.connect(**params)
        cur: psycopg2.extensions.cursor = conn.cursor()

        # dict to track keys with specific database calls
        special_keys: dict[str, int] = {
            "sample_head": 0,
            "sample_file": 0,
            "patient_head": 0,
            "patient_file": 0,
        }
        ref_dir: str = f"{args.ref_dir}/" if args.ref_dir[-1] != "/" else args.ref_dir
        os.makedirs(datasheet_dir, exist_ok=True)
        get_data_clinical(cur, config_data, "sample", ref_dir, datasheet_dir)
        get_data_clinical(cur, config_data, "patient", ref_dir, datasheet_dir)

        # For all other tables to be printed simply, not in special_keys
        for key in config_data["database_pulls"]:
            if key not in special_keys and key != "_comment":
                tbl_name: str = config_data["database_pulls"][key]["table"]
                rows: list[tuple]
                colnames: list[str]
                rows, colnames = generic_pull(cur, tbl_name)
                out_fn: str = config_data["database_pulls"][key]["out_file"]
                # all data_clinical sheets go in this dir
                if out_fn.startswith("data_clinical_"):
                    out_fn = datasheet_dir + "/" + out_fn
                out_file: IO = open(out_fn, "w")
                generic_print(out_file, rows, colnames)
                out_file.close()

        # close the communication with the PostgreSQL
        cur.close()
    except (Exception, psycopg2.DatabaseError) as error:
        print(error, file=sys.stderr)
        conn.close()
    conn.close()


def main():
    parser = argparse.ArgumentParser(
        description="Pull clinical data and genomics file etl support from D3b data warehouse."
    )
    parser.add_argument(
        "-db",
        "--db-ini",
        action="store",
        dest="db_ini",
        help="Database config file - formatting like aws or sbg creds",
    )
    parser.add_argument(
        "-p",
        "--profile",
        action="store",
        dest="profile",
        help="ini profile name",
        default="postgresql",
    )
    parser.add_argument(
        "-sc",
        "--study-config",
        action="store",
        dest="study_config",
        help="json config file with study information; see cbioportal_etl/STUDY_CONFIGS/json_files for examples",
    )
    parser.add_argument(
        "-r",
        "--ref-dir",
        action="store",
        dest="ref_dir",
        help="dir name containing template data_clinical* header files",
    )
    parser.add_argument(
        "-a",
        "--all",
        action="store_true",
        dest="all",
        help="flag to include all relevant files, not just status=active files, NOT RECOMMENDED",
    )

    args = parser.parse_args()
    run_py(args)


if __name__ == "__main__":
    main()
