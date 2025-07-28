"""Check dbt snapshot tables to determine which tables in a specified schema have been updated within a given date range.

Usage:
    python kf-cbioportal-etl/cbioportal_etl/scripts/check_snapshot_updates.py --db-ini /path/to/db.ini

Arguments:
    --schema: The Postgres schema that contains the snapshot tables (default: wongj4_dev_schema_snapshots until prod snapshots are created)
    --db-ini: Path to the db.ini file with Postgres connection details.
    --start:  Start date (inclusive) for detecting updates. Defaults to last 14 days if not provided.
    --end:    End date (inclusive) for detecting updates. Defaults to current time if not provided.

Output:
    - Prints a list of snapshot tables that contain updates during the given date range.
    - Writes the list of updated tables to `updated_snapshots.txt`.

"""

import argparse
import os
import psycopg2

from configparser import ConfigParser
from datetime import datetime, timedelta


def config(filename: str = "database.ini", section: str = "postgresql") -> dict:
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
        e = f"Section {section} not found in the {format} file"
        raise Exception(e)

    return db


def fetch_schema_tables(schema: str, db_ini: str) -> list:
    """Get a list of tables in a schema.
    
    Args: 
        schema (str): Name of schema to query for tables.
        db_ini (str): Path to the database config file.
    Returns:
        tables (list): List of tables in the schema.
    """
    # Load database login info
    params: dict = config(filename=db_ini)
    try:
        conn: psycopg2.extensions.connection = psycopg2.connect(**params)
        cur: psycopg2.extensions.cursor = conn.cursor()
        cur.execute(
            f"""
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = %s
            AND table_type = 'BASE TABLE'
            """, (schema,)
        )
    
        tables = [row[0] for row in cur.fetchall()]
        cur.close()
        conn.close()
        return tables

    except (Exception, psycopg2.DatabaseError) as error:
        print(error)
        return []


def get_updated_tables(schema: str, tables: list, db_ini: str, start_date: datetime, end_date: datetime) -> list:
    """Check which tables have snapshot updates in the date range.

    Args: 
        schema (str): Name of schema to query for tables.
        tables (list): List of table names to check for updates.
        db_ini (str): Path to the database config file.
        start_date (datetime): Start of date range (inclusive).
        end_date (datetime): End of date range (inclusive).
    Returns:
        updated (list): List of tables in the schema with updates in date range.
    """
    params = config(filename=db_ini)
    updated = []

    conn = psycopg2.connect(**params)
    cur = conn.cursor()

    for table in tables:
        query = f"""
            SELECT COUNT(*)
            FROM {schema}.{table}
            WHERE dbt_valid_from BETWEEN %s AND %s
        """
        try:
            cur.execute(query, (start_date, end_date))
            count = cur.fetchone()[0]
            if count > 0:
                updated.append(table)
        except psycopg2.errors.UndefinedColumn:
            conn.rollback()  # reset transaction state
            print(f"Skipping {table}, no dbt_valid_from column")

    cur.close()
    conn.close()
    return updated


def main():
    parser = argparse.ArgumentParser(
        description="Generate dbt snapshot SQL for all *_data_* and *_timeline_* tables in prod_cbio schema."
    )
    parser.add_argument(
        "-db",
        "--db-ini",
        action="store",
        dest="db_ini",
        required=True,
        help="Path to db.ini file with database connection details"
    )
    parser.add_argument(
        "-s",
        "--schema",
        default="wongj4_dev_schema_snapshots",
        required=False,
        help="Postgres schema name to scan."
    )
    parser.add_argument(
        "-sd",
        "--start", 
        required=False,
        help="Start date (YYYY-MM-DD). Default = last 14 days"
    )
    parser.add_argument(
        "-ed",
        "--end", 
        required=False,
        help="End date (YYYY-MM-DD). Default = now"
    )

    args = parser.parse_args()
    now = datetime.now()
    if args.start:
        start_date = datetime.fromisoformat(args.start)
    else:
        start_date = now - timedelta(days=14)

    end_date = datetime.fromisoformat(args.end) if args.end else now

    tables = fetch_schema_tables(args.schema, args.db_ini)
    updated_tables = get_updated_tables(args.schema, tables, args.db_ini, start_date, end_date)

    print(f"Tables with updates between {start_date} and {end_date}:")
    for t in updated_tables:
        print(f"- {t}")

    with open("updated_snapshots.txt", "w") as f:
        for t in updated_tables:
            f.write(t + "\n")


if __name__ == "__main__":
    main()