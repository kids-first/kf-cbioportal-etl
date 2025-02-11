"""
testing DWH SQL

python3 ~/cbio_etl_studies/ETL/generate_config_v2.py \
    --db-ini ~/cbio_etl_studies/ETL/db.ini \
    --study open_chordoma \
    --config-tsv ~/cbio_etl_studies/ETL/test_tsv.tsv

"""

import psycopg2
from psycopg2 import sql
from configparser import ConfigParser
import argparse
import pandas as pd


""" put this in it's own script (to log into database) to import """
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


def get_file_type(db_cur, study, tbl_name):
    """
    Queries 'File_Type' column for dl_file_type_list
    """
    tbl_sql = sql.SQL('SELECT DISTINCT "file_type" FROM prod_cbio.{};').format(sql.Identifier(tbl_name))

    db_cur.execute(tbl_sql)
    file_types = [row[0] for row in db_cur.fetchall()] 

    entries = [] 
    if file_types:
            entries.append([study, "dl_file_type_list", tbl_name, "", ", ".join(file_types)])

    return pd.DataFrame(entries, columns=["Study", "Key", "File", "Sub-Key", "Value"]) if entries else pd.DataFrame(columns=["Study", "Key", "File", "Sub-Key", "Value"])


def get_merged(db_cur, study, tbl_name, config_tsv_df):
    """
    Queries 'File_Type' column, checks for ['cnv', 'rsem', 'fusion', 'maf]
    If exists, create 'merged_' sections 
    """
    tbl_name = f"{study}_genomics_etl_file"
    tbl_sql = sql.SQL('SELECT DISTINCT "File_Type" FROM prod_cbio.{};').format(sql.Identifier(tbl_name))
    
    db_cur.execute(tbl_sql)
    file_types = [row[0].lower() for row in db_cur.fetchall()] 
    
    merged_entries = []
    file_type_mapping = {
        "cnv": "merged_cnvs",
        "maf": "merged_mafs",
        "fusion": "merged_fusion",
        "rsem": "merged_rsem"
    }

    for file_type, merged_key in file_type_mapping.items():
        if file_type in file_types:
            # Add the base 'merged_{type}' row
            merged_entries.append([study, merged_key, f"{study}_case_meta_config.json", "dir", merged_key])
            
            # Extract rows from TSV where Key == "merged_{type}"
            extracted_rows = config_tsv_df[config_tsv_df["Key"] == merged_key].copy()
            extracted_rows["Study"] = study  
            extracted_rows["File"] = f"{study}_case_meta_config.json"  
            
            merged_entries.extend(extracted_rows.values.tolist())

    return pd.DataFrame(merged_entries, columns=["Study", "Key", "File", "Sub-Key", "Value"]) if merged_entries else pd.DataFrame(columns=["Study", "Key", "File", "Sub-Key", "Value"])


def get_dna_rna_ext_list(db_cur, study, tbl_name):
    tbl_sql = sql.SQL("""SELECT DISTINCT substring("File_Name" FROM position ('.' IN "File_Name") + 1) AS file_name_ext, "file_type" FROM prod_cbio.{};""").format(sql.Identifier(tbl_name))
    
    db_cur.execute(tbl_sql)
    file_exts = [row[0].lower() for row in db_cur.fetchall()]

    entries = []
    
    file_type_mapping = {
        "rna_ext_list": {
            "expression": "rsem.genes.results.gz",
            "fusion": "annofuse_filter.tsv"
        },
        "dna_ext_list": {
            "seg": "controlfreec.seg"
        }
    }

    # RNA & DNA file checks
    for key, sub_keys in file_type_mapping.items():
        for sub_key, file_name in sub_keys.items():
            if file_name in file_exts:
                entries.append([study, key, f"{study}_data_processing_config.json", sub_key, file_name])

    # Check for "CNVs.p.value.txt" as a suffix (not an exact match)
    has_cnv = any(ext.endswith("cnvs.p.value.txt") for ext in file_exts)
    if has_cnv:
        entries.append([study, "dna_ext_list", f"{study}_data_processing_config.json", "copy_number", "CNVs.p.value.txt"])

    # Handle mutation separately since there may be multiple `.maf` files depending on the study
    maf_files = sorted({ext for ext in file_exts if ext.endswith(".maf")})
    if maf_files:
        entries.append([study, "dna_ext_list", f"{study}_data_processing_config.json", "mutation", ", ".join(maf_files)])

    return pd.DataFrame(entries, columns=["Study", "Key", "File", "Sub-Key", "Value"]) if entries else pd.DataFrame(columns=["Study", "Key", "File", "Sub-Key", "Value"])


def get_file_loc_defs(db_cur, study, tbl_name):
    """
    Queries 'file_type' column for file_loc_defs
    """
    tbl_sql = sql.SQL("""SELECT DISTINCT substring("File_Name" FROM position ('.' IN "File_Name") + 1) AS file_name_ext, "file_type" FROM prod_cbio.{};""").format(sql.Identifier(tbl_name))

    db_cur.execute(tbl_sql)
    rows = db_cur.fetchall()
    df = pd.DataFrame(rows, columns=["file_name_ext", "file_type"])

    entries = []

    mappings = {
        "ctrlfreec_info": ("cnvs.info", "ctrlfreec_info"),
        "ctrlfreec_pval": ("cnvs.pval", "ctrlfreec_pval"),
        "ctrlfreec_bam_seg": ("cnvs.seg", "ctrlfreec_bam_seg"),
        "DGD_FUSION": ("dgd_fusion", "fusion-dgd.tsv.gz"),  
        "annofuse_filtered_fusions_tsv": ("fusion", "annofuse_filtered_fusions_tsv"),
        "RSEM_gene": ("rsem", "RSEM_gene")
    }

    for file_type, (sub_key, output_value) in mappings.items():
        if file_type in df["file_type"].values:
            entries.append([study, "file_loc_defs", f"{study}_data_processing_config.json", sub_key, output_value])

    # use file_name_ext to get mafs
    maf_files = df[df["file_name_ext"].str.endswith(".maf", na=False)]["file_type"].unique().tolist()
    if maf_files:
        entries.append([study, "file_loc_defs", f"{study}_data_processing_config.json", "mafs.kf", ", ".join(maf_files)])

    return pd.DataFrame(entries, columns=["Study", "Key", "File", "Sub-Key", "Value"]) if entries else pd.DataFrame(columns=["Study", "Key", "File", "Sub-Key", "Value"])


def run_py(args):
    # Load database login info
    params = config(filename=args.db_ini, section=args.profile)
    study = args.study
    tbl_name = f"{study}_genomics_etl_file"

    try:
        conn = psycopg2.connect(**params)
        cur = conn.cursor()

        config_tsv_df = pd.read_csv(args.config_tsv, sep="\t")

        file_type_df = get_file_type(cur, study, tbl_name)
        merged_types_df = get_merged(cur, study, tbl_name, config_tsv_df)
        dna_ext_list_df = get_dna_rna_ext_list(cur, study, tbl_name)
        file_loc_defs_df = get_file_loc_defs(cur, study, tbl_name)

        final_df = pd.concat([file_type_df, merged_types_df, dna_ext_list_df, file_loc_defs_df], ignore_index=True)
        print(final_df)

        cur.close()

    except (Exception, psycopg2.DatabaseError) as error:
        print(error)
        conn.close()
    conn.close()


def main():
    parser = argparse.ArgumentParser(description="Pull clinical data and genomics file etl support from D3b data warehouse.")
    parser.add_argument("-db", "--db-ini", action="store", dest="db_ini", help="Database config file - formatting like aws or sbg creds")
    parser.add_argument("-p", "--profile", action="store", dest="profile", help="ini profile name", default="postgresql")
    parser.add_argument("-s", "--study", action="store", dest="study", help="Cancer study ID to compare on server")
    parser.add_argument("-v", "--config-tsv", action="store", dest="config_tsv", help="Path to the input TSV file")

    args = parser.parse_args()
    run_py(args)

if __name__ == "__main__":
    main()