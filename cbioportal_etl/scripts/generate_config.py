#!/usr/bin/env python3
"""Queries DWH and outputs values to TSV.

Generate a JSON config from a TSV based on the provided study.

Command:
python3 /kf-cbioportal-etl/cbioportal_etl/cbioportal_etl/scripts/generate_config.py \
    --db-ini db.ini \
    --study study \
    --config-tsv /kf-cbioportal-etl/cbioportal_etl/STUDY_CONFIGS/all_studies_config_values.tsv
"""

import argparse
import ast
import json
from configparser import ConfigParser

import pandas as pd
import psycopg2
from psycopg2 import sql


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


def get_file_type(
    db_cur: psycopg2.extensions.cursor, study: str, tbl_name: str, schema: str) -> pd.DataFrame:
    """Get valid list of 'file_type' from genomics_etl_file table.

    Args:
        db_cur (psycopg2.extensions.cursor): Database cursor
        study (str): Study ID
        tbl_name (str): Table name to query
        schema (str): Schema name, default is 'prod_cbio'

    Returns:
        pd.DataFrame: DataFrame with file types

    """
    tbl_sql = sql.SQL('SELECT DISTINCT "file_type" FROM {}.{};').format(
        sql.Identifier(schema) , sql.Identifier(tbl_name)
    )

    db_cur.execute(tbl_sql)
    file_types = [row[0] for row in db_cur.fetchall()]

    entries = []
    if file_types:
        entries.append([study, "dl_file_type_list", tbl_name, "", ", ".join(file_types)])

    return (
        pd.DataFrame(entries, columns=["Study", "Key", "File", "Sub-Key", "Value"])
        if entries
        else pd.DataFrame(columns=["Study", "Key", "File", "Sub-Key", "Value"])
    )


def get_merged(db_cur: psycopg2.extensions.cursor, study: str, tbl_name: str, config_tsv_df: pd.DataFrame, schema: str) -> pd.DataFrame:
    """Query 'etl_file_type' column, checks for ['cnv', 'rsem', 'fusion', 'maf].

    If exists, create 'merged_' sections

    Args:
        db_cur (psycopg2.extensions.cursor): Database cursor
        study (str): Study ID
        tbl_name (str): Table name to query
        config_tsv_df (pd.DataFrame): DataFrame containing static config values
        schema (str): Schema name, default is 'prod_cbio'

    Returns:
        pd.DataFrame: DataFrame with merged entries

    """
    tbl_name = f"{study}_genomics_etl_file"
    tbl_sql = sql.SQL('SELECT DISTINCT "etl_file_type" FROM {}.{};').format(
        sql.Identifier(schema) , sql.Identifier(tbl_name)
    )

    db_cur.execute(tbl_sql)
    file_types = [row[0].lower() for row in db_cur.fetchall()]

    merged_entries = []
    file_type_mapping = {
        "cnv": "merged_cnvs",
        "maf": "merged_mafs",
        "fusion": "merged_fusion",
        "rsem": "merged_rsem",
    }

    for file_type, merged_key in file_type_mapping.items():
        if file_type in file_types:
            # Add the base 'merged_{type}' row
            merged_entries.append(
                [study, merged_key, f"{study}_case_meta_config.json", "dir", merged_key]
            )

            # Extract rows from TSV where Key == "merged_{type}"
            extracted_rows = config_tsv_df[
                ((config_tsv_df["Key"] == merged_key) & (config_tsv_df["Study"] == study)) |
                ((config_tsv_df["Key"] == merged_key) & (config_tsv_df["Study"] == f"{file_type}_etl_file_type_only"))
            ].copy()

            merged_entries.extend(extracted_rows.values.tolist())

    return (
        pd.DataFrame(merged_entries, columns=["Study", "Key", "File", "Sub-Key", "Value"])
        if merged_entries
        else pd.DataFrame(columns=["Study", "Key", "File", "Sub-Key", "Value"])
    )


def get_dna_rna_ext_list(db_cur: psycopg2.extensions.cursor, study: str, tbl_name: str, schema: str) -> pd.DataFrame:
    """Query 'file_name' column for dna_ext_list and rna_ext_list.

    Args:
        db_cur (psycopg2.extensions.cursor): Database cursor
        study (str): Study ID
        tbl_name (str): Table name to query
        schema (str): Schema name, default is 'prod_cbio'
    Returns:
        pd.DataFrame: DataFrame with file extensions and types

    """
    tbl_sql = sql.SQL(
        """SELECT DISTINCT substring("file_name" FROM position ('.' IN "file_name") + 1) AS file_name_ext, "file_type" FROM {}.{};"""
    ).format(sql.Identifier(schema) , sql.Identifier(tbl_name))

    db_cur.execute(tbl_sql)
    file_exts = [row[0] for row in db_cur.fetchall()]

    entries = []

    file_type_mapping = {
        "rna_ext_list": {"expression": "rsem.genes.results.gz", "fusion": "annoFuse_filter.tsv"},
        "dna_ext_list": {"seg": "controlfreec.seg"},
    }

    # RNA & DNA file checks
    for key, sub_keys in file_type_mapping.items():
        for sub_key, file_name in sub_keys.items():
            if file_name in file_exts:
                entries.append(
                    [study, key, f"{study}_data_processing_config.json", sub_key, file_name]
                )

    # Check for "CNVs.p.value.txt" as a suffix (not an exact match)
    has_cnv = any(ext.endswith("CNVs.p.value.txt") for ext in file_exts)
    if has_cnv:
        entries.append(
            [
                study,
                "dna_ext_list",
                f"{study}_data_processing_config.json",
                "copy_number",
                "CNVs.p.value.txt",
            ]
        )

    # Handle mutation separately since there may be multiple `.maf` files depending on the study
    maf_files = sorted({ext for ext in file_exts if ext.endswith(".maf")})
    if maf_files:
        entries.append(
            [
                study,
                "dna_ext_list",
                f"{study}_data_processing_config.json",
                "mutation",
                ", ".join(maf_files),
            ]
        )

    return (
        pd.DataFrame(entries, columns=["Study", "Key", "File", "Sub-Key", "Value"])
        if entries
        else pd.DataFrame(columns=["Study", "Key", "File", "Sub-Key", "Value"])
    )


def get_file_loc_defs(db_cur: psycopg2.extensions.cursor, study: str, tbl_name: str, schema: str) -> pd.DataFrame:
    """Get file name extension 'file_type' column for file_loc_defs.

    Map file extensions to file_type to inform proper grouping of files by file_type

    Args:
        db_cur (psycopg2.extensions.cursor): Database cursor
        study (str): Study ID
        tbl_name (str): Table name to query
        schema (str): Schema name, default is 'prod_cbio'
    Returns:
        pd.DataFrame: DataFrame with file location definitions
    """
    tbl_sql = sql.SQL(
        """SELECT DISTINCT substring("file_name" FROM position ('.' IN "file_name") + 1) AS file_name_ext, "file_type" FROM {}.{};"""
    ).format(sql.Identifier(schema) , sql.Identifier(tbl_name))

    db_cur.execute(tbl_sql)
    rows = db_cur.fetchall()
    df = pd.DataFrame(rows, columns=["file_name_ext", "file_type"])

    entries = []

    mappings = {
        "ctrlfreec_info": ("cnvs.info", "ctrlfreec_info"),
        "ctrlfreec_pval": ("cnvs.pval", "ctrlfreec_pval"),
        "ctrlfreec_bam_seg": ("cnvs.seg", "ctrlfreec_bam_seg"),
        "annofuse_filtered_fusions_tsv": ("fusion", "annofuse_filtered_fusions_tsv"),
        "RSEM_gene": ("rsem", "RSEM_gene"),
    }

    for file_type, (sub_key, output_value) in mappings.items():
        if file_type in df["file_type"].values:
            entries.append(
                [
                    study,
                    "file_loc_defs",
                    f"{study}_data_processing_config.json",
                    sub_key,
                    output_value,
                ]
            )

    # use file_name_ext to get mafs
    maf_files = (
        df[df["file_name_ext"].str.endswith(".maf", na=False)]["file_type"].unique().tolist()
    )
    if maf_files:
        entries.append(
            [
                study,
                "file_loc_defs",
                f"{study}_data_processing_config.json",
                "mafs.kf",
                ", ".join(maf_files),
            ]
        )

    return (
        pd.DataFrame(entries, columns=["Study", "Key", "File", "Sub-Key", "Value"])
        if entries
        else pd.DataFrame(columns=["Study", "Key", "File", "Sub-Key", "Value"])
    )


def get_cases_info(db_cur: psycopg2.extensions.cursor, study: str, tbl_name: str, schema: str) -> pd.DataFrame:
    """Query DWH for study-specific case lists.

    - 'cases_all' (all tumor or tumor & model samples)
    - 'cases_3way_complete' (only if the study has mutation, CNA, and mRNA data)
    - 'cases_cna' (copy number data)
    - 'cases_cnaseq' (mutation and CNA data)
    """
    # List of studies with both Tumor and Model samples
    tumor_model_studies = {
        "bllnos_sd_z6mwd3h0_2018",
        "chdm_phs001643_2018",
        "chdm_phs002301_2021",
        "chdm_sd_7spqtt8m",
        "open_chordoma",
        "os_sd_zxjffmef_2015",
    }
    is_tumor_model = study in tumor_model_studies
    original_config = f"{study}_case_meta_config.json"

    # Set case_list descriptions
    case_list_desc = "All tumor and model samples" if is_tumor_model else "All tumor samples"
    case_list_name = "All tumors and models" if is_tumor_model else "All tumors"
    # cases_3way_complete
    case_list_desc_3way = f"{case_list_desc} with mutation, CNA, and mRNA data"
    case_list_name_3way = (
        "Tumor and model samples with mutation, CNA and mRNA data"
        if is_tumor_model
        else "Tumor samples with mutation, CNA and mRNA data"
    )
    # cases_cna
    case_list_desc_cna = f"{case_list_desc} with CNA data"
    case_list_name_cna = (
        "Tumor and model samples with CNA data" if is_tumor_model else "Tumor Samples with CNA data"
    )
    # cases_cnaseq
    case_list_desc_cnaseq = "All tumor samples with mutation and CNA data"
    case_list_name_cnaseq = (
        "Tumor samples with mutation and CNA data"
        if is_tumor_model
        else "Tumor samples with mutation and CNA data"
    )
    # cases_RNA_Seq_v2_mRNA
    case_list_desc_rna = "All samples with mRNA expression data"
    case_list_name_rna = (
        "Tumor and model samples with mRNA data (RNA Seq V2)"
        if is_tumor_model
        else "Tumor Samples with mRNA data (RNA Seq V2)"
    )
    # cases_sequenced
    case_list_desc_sequenced = f"{case_list_desc} with mutation data"
    case_list_name_sequenced = (
        "Tumor and model samples with mutations"
        if is_tumor_model
        else "Tumor samples with mutations"
    )
    # cases_sv
    case_list_desc_sv = f"{case_list_desc} with fusion data"
    case_list_name_sv = (
        "Tumor and model samples with fusions" if is_tumor_model else "Tumor samples with fusions"
    )

    # Create 'cases_all'
    entries = [
        [study, "cases_all", original_config, "case_list_description", case_list_desc],
        [study, "cases_all", original_config, "case_list_name", case_list_name],
    ]

    # Query for other cases
    tbl_sql = sql.SQL(
        """SELECT DISTINCT substring("file_name" FROM position ('.' IN "file_name") + 1) AS file_name_ext, "file_type" FROM {}.{};"""
    ).format(sql.Identifier(schema) , sql.Identifier(tbl_name))
    db_cur.execute(tbl_sql)
    rows = db_cur.fetchall()
    file_name_exts = {row[0] for row in rows if row[0] is not None}

    has_maf = any(name.endswith(".maf") for name in file_name_exts)
    has_cnv = any("CNVs.p.value.txt" in name for name in file_name_exts)
    has_rsem = "rsem.genes.results.gz" in file_name_exts
    has_sv = "annoFuse_filter.tsv" in file_name_exts

    # cases_3way_complete - check if study has maf, cnv, and rsem
    if has_maf and has_cnv and has_rsem:
        entries.extend([
            [study, "cases_3way_complete", original_config, "case_list_category", "all_cases_with_mutation_and_cna_and_mrna_data"],
            [study, "cases_3way_complete", original_config, "case_list_description", case_list_desc_3way],
            [study, "cases_3way_complete", original_config, "case_list_name", case_list_name_3way],
            [study, "cases_3way_complete", original_config, "stable_id", "3way_complete"]
        ])
    # cna - check if study has cnv
    if has_cnv:
        entries.extend([
            [study, "cases_cna", original_config, "case_list_category", "all_cases_with_cna_data"],
            [study, "cases_cna", original_config, "case_list_description", case_list_desc_cna],
            [study, "cases_cna", original_config, "case_list_name", case_list_name_cna],
            [study, "cases_cna", original_config, "stable_id", "cna"]
        ])
    # cnaseq - check if study has maf + cnv
    if has_maf and has_cnv:
        entries.extend([
            [study, "cases_cnaseq", original_config, "case_list_category", "all_cases_with_mutation_and_cna_data"],
            [study, "cases_cnaseq", original_config, "case_list_description", case_list_desc_cnaseq],
            [study, "cases_cnaseq", original_config, "case_list_name", case_list_name_cnaseq],
            [study, "cases_cnaseq", original_config, "stable_id", "cnaseq"]
        ])
    # cases_RNA_Seq_v2_mRNA - check if study has rsem
    if has_rsem:
        entries.extend([
            [study, "cases_RNA_Seq_v2_mRNA", original_config, "case_list_category", "all_cases_with_mrna_rnaseq_data"],
            [study, "cases_RNA_Seq_v2_mRNA", original_config, "case_list_description", case_list_desc_rna],
            [study, "cases_RNA_Seq_v2_mRNA", original_config, "case_list_name", case_list_name_rna],
            [study, "cases_RNA_Seq_v2_mRNA", original_config, "stable_id", "rna_seq_v2_mrna"]
        ])
    # cases_sequenced - check if study has maf
    if has_maf:
        entries.extend([
            [study, "cases_sequenced", original_config, "case_list_category", "all_cases_with_mutation_data"],
            [study, "cases_sequenced", original_config, "case_list_description", case_list_desc_sequenced],
            [study, "cases_sequenced", original_config, "case_list_name", case_list_name_sequenced],
            [study, "cases_sequenced", original_config, "stable_id", "sequenced"]
        ])
    # cases_sv - check if study has fusion data
    if has_sv:
        entries.extend([
            [study, "cases_sv", original_config, "case_list_category", "all_cases_with_sv_data"],
            [study, "cases_sv", original_config, "case_list_description", case_list_desc_sv],
            [study, "cases_sv", original_config, "case_list_name", case_list_name_sv],
            [study, "cases_sv", original_config, "stable_id", "sv"]
        ])

    return pd.DataFrame(entries, columns=["Study", "Key", "File", "Sub-Key", "Value"])


def get_database_pulls(db_cur: psycopg2.extensions.cursor, study: str, tbl_name: str, schema: str) -> pd.DataFrame:
    """Query DWH for study-specific timeline tables.

    If a study has *_timeline_clinical_event, *_timeline_imaging, *_timeline_specimen, *_timeline_surgery,
    or *_timeline_treatment tables, add database_pulls and data_sheets entries.

    Also, adds:
    - genomics_etl.table for all studies
    - patient_file.table for all studies
    - sample_file.table for all studies

    Args:
        db_cur (psycopg2.extensions.cursor): Database cursor
        study (str): Study ID
        tbl_name (str): Table name to query
        schema (str): Schema name, default is 'prod_cbio'
    Returns:
        pd.DataFrame: DataFrame with database pulls and data sheets entries
    """
    original_config = f"{study}_case_meta_config.json"
    entries = ([
        [study, "database_pulls", original_config, "genomics_etl.table", f"prod_cbio.{tbl_name}"],
        [study, "database_pulls", original_config, "patient_file.table", f"prod_cbio.{study}_data_clinical_patient"],
        [study, "database_pulls", original_config, "sample_file.table", f"prod_cbio.{study}_data_clinical_sample"]   
    ])

    # Query to check if the timeline tables exist for the study
    timeline_tables = {
        "clinical_event_timeline": {
            "table": f"{study}_timeline_clinical_event",
            "out_file": "data_clinical_timeline_clinical_event.txt",
        },
        "imaging_timeline": {
            "table": f"{study}_timeline_imaging",
            "out_file": "data_clinical_timeline_imaging.txt",
        },
        "specimen_timeline": {
            "table": f"{study}_timeline_specimen",
            "out_file": "data_clinical_timeline_specimen.txt",
        },
        "surgery_timeline": {
            "table": f"{study}_timeline_surgery",
            "out_file": "data_clinical_timeline_surgery.txt",
        },
        "treatment_timeline": {
            "table": f"{study}_timeline_treatment",
            "out_file": "data_clinical_timeline_treatment.txt",
        },
    }

    timeline_subkey_map = {
        "clinical_event_timeline": "dtypes.clinical_event",
        "imaging_timeline": "dtypes.imaging",
        "specimen_timeline": "dtypes.specimen",
        "surgery_timeline": "dtypes.surgery",
        "treatment_timeline": "dtypes.treatment",
    }

    for key, info in timeline_tables.items():
        check_table_sql = sql.SQL("SELECT to_regclass('{}.{}') IS NOT NULL;").format(
            sql.Identifier(schema),
            sql.Identifier(info["table"])
        )
        db_cur.execute(check_table_sql)
        table_exists = db_cur.fetchone()[0]

        if table_exists:
            entries.extend([
                [study, "database_pulls", original_config, f"{key}.out_file", info["out_file"]],
                [study, "database_pulls", original_config, f"{key}.table", f"prod_cbio.{info['table']}"]
            ])

            correct_subkey = timeline_subkey_map.get(key, f"dtypes.{key}")

            entries.append(
                [
                    study,
                    "data_sheets",
                    original_config,
                    correct_subkey,
                    json.dumps(
                        {
                            "_comment": "see https://docs.cbioportal.org/file-formats/#event-types for detailed specifics",
                            "cbio_name": info["out_file"],
                            "meta_file_attr": {
                                "genetic_alteration_type": "CLINICAL",
                                "datatype": "TIMELINE",
                            },
                        }
                    ),
                ]
            )

    return pd.DataFrame(entries, columns=["Study", "Key", "File", "Sub-Key", "Value"])


def parse_value(value: str, should_split: bool=False, sub_key: str | None=None) -> str | int | list | dict:
    """Parse the value based on its type.

    Args:
        value (str): The value to parse.
        should_split (bool): Whether to split the value into a list.
        sub_key (str | None): The sub-key to determine specific parsing behavior.

    Returns:
        The parsed value, which can be a string, int, list, or dict.
    """
    value = value.strip()

    if sub_key == "pmid":
        return value

    if should_split:
        split_values = [v.strip() for v in value.split(",")]
        return split_values if len(split_values) > 1 else split_values[0]
    try:
        return json.loads(value.replace("'", '"'))
    except json.JSONDecodeError:
        try:
            return ast.literal_eval(value)
        except (ValueError, SyntaxError):
            return int(value) if value.isdigit() else value


def update_nested_dict(result: dict, key: str, sub_key: str, value: str) -> None:
    """Update a nested dictionary with the given key, sub-key, and value.

    Args:
        result (dict): The dictionary to update.
        key (str): The main key for the dictionary.
        sub_key (str): The sub-key for the nested dictionary.
        value (str): The value to set for the key or sub-key.

    """
    should_split = (sub_key and sub_key.lower() in {"mutation", "mafs.kf"}) or (
        key and key.lower() in {"dl_file_type_list"}
    )
    current = result.setdefault(key, {})
    if sub_key:
        keys = sub_key.split(".")
        for part in keys[:-1]:
            current = current.setdefault(part, {})
        current[keys[-1]] = parse_value(value, should_split, sub_key)
    else:
        result[key] = parse_value(value, should_split, sub_key=None)


def generate_json(study_df: pd.DataFrame, study: str):
    """Generate JSON file from study DataFrame.

    Args:
        study_df (pd.DataFrame): DataFrame containing study configuration values.
        study (str): Study ID for naming the output JSON file.
    """
    result = {}

    for _, row in study_df.iterrows():
        key = row["Key"]
        sub_key = row["Sub-Key"] if not pd.isna(row["Sub-Key"]) else ""
        value = row["Value"]

        update_nested_dict(result, key, sub_key, value)

    output_file = f"{study}_config.json"
    with open(output_file, "w") as json_file:
        json.dump(result, json_file, indent=4)

    print(f"JSON file generated: {output_file}")


def generate_tsv(args: argparse.Namespace) -> None:
    """Generate study-specific TSV.

    Uses a source TSV to created a study-specific TSV

    """
    # Load database login info
    params: dict = config(filename=args.db_ini, section=args.profile)
    study: str = args.study
    schema: str = args.schema
    tbl_name: str = f"{study}_genomics_etl_file"

    studies_without_data_configs: list[str] = [
        "case_cptac",
        "chdm_phs002301_2021",
        "dgd_non_brain",
        "tll_sd_aq9kvn5p_2019",
        "x01_nbl_16_maris",
    ]

    try:
        conn: psycopg2.extensions.connection = psycopg2.connect(**params)
        cur: psycopg2.extensions.cursor = conn.cursor()

        config_tsv_df: pd.DataFrame = pd.read_csv(args.config_tsv, sep="\t")
        static_df: pd.DataFrame = config_tsv_df[
            (config_tsv_df["Study"] == "static")
            | ((config_tsv_df["Study"] == study) & ~config_tsv_df["Key"].str.startswith("merged_"))
        ]

        file_type_df = get_file_type(cur, study, tbl_name, schema)
        merged_types_df = get_merged(cur, study, tbl_name, config_tsv_df, schema)
        dna_ext_list_df = get_dna_rna_ext_list(cur, study, tbl_name, schema)
        file_loc_defs_df = get_file_loc_defs(cur, study, tbl_name, schema)
        cases_df = get_cases_info(cur, study, tbl_name, schema)
        database_pulls_df = get_database_pulls(cur, study, tbl_name, schema)

        final_df = pd.concat(
            [
                static_df,
                file_type_df,
                merged_types_df,
                dna_ext_list_df,
                file_loc_defs_df,
                cases_df,
                database_pulls_df,
            ],
            ignore_index=True,
        )
        if study in studies_without_data_configs:
            extra_entries: list[list[str]] = [
                ["data_processing_config", "cna_flag", "data_processing_config.json", "", "1"],
                ["data_processing_config", "rna_flag", "data_processing_config.json", "", "1"],
            ]
            extra_df: pd.DataFrame = pd.DataFrame(
                extra_entries, columns=["Study", "Key", "File", "Sub-Key", "Value"]
            )
            final_df: pd.DataFrame = pd.concat([final_df, extra_df], ignore_index=True)
        final_df = final_df.sort_values(
            by=["Study", "Key", "Sub-Key"], ascending=[False, True, True]
        )
        final_df.to_csv(f"{study}_config_values.tsv", sep="\t", index=False)
        print(f"TSV file generated: {study}_config_values.tsv.")

        generate_json(final_df, study)

        cur.close()

    except (Exception, psycopg2.DatabaseError) as error:
        print(error)
        conn.close()
    conn.close()


def run_py(args: argparse.Namespace) -> None:
    """If a user provided an existing study TSV, skip TSV generation and just create JSON.

    Otherwise, generate the TSV from scratch
    """
    if args.study_tsv:
        print(f"Loading provided TSV: {args.study_tsv} and converting to JSON...")
        study_tsv = pd.read_csv(args.study_tsv, sep="\t")
        generate_json(study_tsv, args.study)
    else:
        generate_tsv(args)


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
        "-s", "--study", action="store", dest="study", help="Cancer study ID to compare on server"
    )
    parser.add_argument(
        "-ctsv",
        "--config-tsv",
        action="store",
        dest="config_tsv",
        help="Path to the input TSV file containing static values (preset in repo)",
    )
    parser.add_argument(
        "-stsv",
        "--study-tsv",
        action="store",
        dest="study_tsv",
        help="Path to the input TSV file containing completed values for json generation",
    )
    parser.add_argument(
        "--schema",
        action="store",
        dest="schema",
        help="Set to change default schema for table searches for dev purposes only",
        default="prod_cbio",
    )

    args = parser.parse_args()
    run_py(args)


if __name__ == "__main__":
    main()
