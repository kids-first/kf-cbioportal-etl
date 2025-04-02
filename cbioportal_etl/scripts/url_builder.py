#!/usr/bin/env python3
"""
Generates a PedcBioPortal URL with embedded JSON filters based on user-selected study and filters.
This script can be used interactively or via CLI arguments.
Usage: Optimized for pbta_all
"""
import argparse
import json
import psycopg2
from psycopg2 import sql
from configparser import ConfigParser

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


def get_subcohort_selection(db_cur, study):
    """
    Queries {study}_data_clinical_sample for unique values in SUB_COHORT column
    """
    table_name = f"{study}_data_clinical_sample"
    tbl_sql = sql.SQL('SELECT DISTINCT "SUB_COHORT" FROM prod_cbio.{};').format(sql.Identifier(table_name))

    db_cur.execute(tbl_sql)
    subcohort_options = [row[0] for row in db_cur.fetchall()] 
    print("\nAvailable SUB_COHORT(s):")
    for idx, name in enumerate(subcohort_options, start=1):
        print(f"{idx}. {name}")
    print("0. All")

    selected = input("Enter the number(s) of the SUB_COHORT(s) to include (e.g. 1,3 or 0 for All): ").strip()

    if selected == "0":
        return subcohort_options

    try:
        indices = [int(i) for i in selected.split(",")]
        selected_subcohorts = [subcohort_options[i - 1] for i in indices if 1 <= i <= len(subcohort_options)]
        return selected_subcohorts
    except Exception as e:
        print(f"Invalid input: {e}")
        return get_subcohort_selection(db_cur, study)


def prompt_attribute_selection(all_attributes):
    """
    Provides list of attributes user can choose from for filtering
    """
    print("\nAvailable filterable attributes:")
    for i, attribute in enumerate(all_attributes, start=1):
        print(f"{i}. {attribute}")
    print("0. All")
    print("-1. None")

    selected = input("Select attributes to filter by (e.g. 1,3,5 or 0 for All or -1 for None): ").strip()

    if selected == "0":
        return all_attributes
    if selected == "-1":
        return []

    try:
        indices = [int(i) for i in selected.split(",")]
        selected_attrs = [all_attributes[i - 1] for i in indices if 1 <= i <= len(all_attributes)]
        return selected_attrs
    except Exception as e:
        print(f"Invalid input: {e}")
        return prompt_attribute_selection(all_attributes)


def get_attribute_values(cur, study, sample_attributes, selected_attributes, selected_subcohorts):
    """
    Queries {study}_data_clinical_sample and {study}_data_clinical_patient tables
    Lists available attribute values from patient and sample tables
    Only returns values for samples in the selected subcohorts
    """
    attribute_options = {}
    sample_table = f"{study}_data_clinical_sample"
    patient_table = f"{study}_data_clinical_patient"

    for attribute in selected_attributes:
        table = "sample" if attribute in sample_attributes else "patient"
        if table == "sample":
            cur.execute(
                sql.SQL(f"""
                    SELECT DISTINCT "{attribute}" FROM prod_cbio.{sample_table}
                    WHERE "SUB_COHORT" IN %s AND "{attribute}" IS NOT NULL
                    ORDER BY "{attribute}";
                """), (tuple(selected_subcohorts),)
            )
        else:
            cur.execute(
                sql.SQL(f"""
                    SELECT DISTINCT p."{attribute}" FROM prod_cbio.{patient_table} p
                    JOIN prod_cbio.{sample_table} s ON p."PATIENT_ID" = s."PATIENT_ID"
                    WHERE s."SUB_COHORT" IN %s AND p."{attribute}" IS NOT NULL
                    ORDER BY p."{attribute}";
                """), (tuple(selected_subcohorts),)
            )
        values = [row[0] for row in cur.fetchall()]
        if values:
            attribute_options[attribute] = values

    return attribute_options


def prompt_attribute_values(attribute_options):
    """
    Choose values for each attribute to filter for
    """
    selected_filters = []

    for attribute, values in attribute_options.items():
        print(f"\nAvailable values for {attribute}:")
        for i, val in enumerate(values, start=1):
            print(f"{i}. {val}")
        print("0. All")

        selected = input(f"Select values for {attribute} (e.g. 1,2 or 0 for All): ").strip()

        if selected == "0":
            selected_filters.append({"attributeId": attribute, "values": [{"value": val} for val in values]})
        else:
            try:
                indices = [int(i) for i in selected.split(",")]
                selected_values = [{"value": values[i - 1]} for i in indices if 1 <= i <= len(values)]
                if selected_values:
                    selected_filters.append({"attributeId": attribute, "values": selected_values})
            except Exception as e:
                print(f"Invalid input: {e}. Skipping {attribute}.")

    return selected_filters


def prompt_gene_filters(study):
    """
    Prompts user to input a comma-separated list of genes if desired
    """
    gene_types = {
        "Mutated Genes": f"{study}_mutations",
        "CNA Genes": f"{study}_cna",
        "Structural Variant Genes": f"{study}_structural_variants"
    }

    options = list(gene_types.keys())

    print("\nWhich gene categories would you like to filter by?")
    for i, option in enumerate(options, start=1):
        print(f"{i}. {option}")
    print("0. None")

    selected = input("Select gene types (e.g. 1,3 or 0 for None): ").strip()

    if selected == "0":
        return []

    try:
        indices = [int(i) for i in selected.split(",")]
        selected_types = [options[i - 1] for i in indices if 1 <= i <= len(options)]
    except Exception as e:
        print(f"Invalid input: {e}")
        return prompt_gene_filters(study)

    gene_filters = []
    for gene_type in selected_types:
        gene_input = input(f"Enter comma-separated gene symbols for {gene_type}: ").strip()
        gene_list = [g.strip() for g in gene_input.split(",") if g.strip()]
        if not gene_list:
            continue

        gene_queries = [[{
            "hugoGeneSymbol": gene,
            "entrezGeneId": 0,
            "alterations": ["AMP", "HOMDEL"] if "CNA" in gene_type else [],
            "includeDriver": True,
            "includeVUS": True,
            "includeUnknownOncogenicity": True,
            "tiersBooleanMap": {},
            "includeUnknownTier": True,
            "includeGermline": True,
            "includeSomatic": True,
            "includeUnknownStatus": True
        } for gene in gene_list]]

        gene_filters.append({
            "molecularProfileIds": [gene_types[gene_type]],
            "geneQueries": gene_queries
        })

    return gene_filters


def build_json(study, sub_cohorts, clinical_filters, gene_filters):
    """
    Uses all user-selected filters to create a JSON that will be embedded in URL
    """
    filter_json = {
        "studyIds": [study],
        "clinicalDataFilters": [{"attributeId": "SUB_COHORT", "values": [{"value": val} for val in sub_cohorts]}] + clinical_filters,
        "alterationFilter": {
            "copyNumberAlterationEventTypes": {"AMP": True, "HOMDEL": True},
            "mutationEventTypes": {"any": True},
            "structuralVariants": None,
            "includeDriver": True,
            "includeVUS": True,
            "includeUnknownOncogenicity": True,
            "includeUnknownTier": True,
            "includeGermline": True,
            "includeSomatic": True,
            "includeUnknownStatus": True,
            "tiersBooleanMap": {}
        }
    }
    if gene_filters:
        filter_json["geneFilters"] = gene_filters

    return filter_json


def run_py(args):
    params = config(filename=args.db_ini, section=args.profile)
    study = args.study
    try:
        conn = psycopg2.connect(**params)
        cur = conn.cursor()
        # Choose SUB_COHORT
        sub_cohorts = get_subcohort_selection(cur, study)
        # Choose attributes to filter and values to filter for
        patient_attributes = ["SEX"]
        sample_attributes = [
            "SPECIMEN_ID", "EXPERIMENT_STRATEGY", "CANCER_TYPE", "CANCER_TYPE_DETAILED",
            "ONCOTREE_CODE", "TUMOR_TISSUE_TYPE", "BROAD_HISTOLOGY", "CBTN_TUMOR_TYPE",
            "TUMOR_TYPE", "SAMPLE_TYPE", "CANCER_GROUP", "RNA_LIBRARY_SELECTION"
            ]
        all_attributes = patient_attributes + sample_attributes
        selected_attributes = prompt_attribute_selection(all_attributes)
        attribute_options = get_attribute_values(cur, study, sample_attributes, selected_attributes, sub_cohorts)
        clinical_filters = prompt_attribute_values(attribute_options)
        gene_filters = prompt_gene_filters(study)
        # Generate URL
        filter_json = build_json(study, sub_cohorts, clinical_filters, gene_filters)
        encoded_filter = json.dumps(filter_json)
        url = f"https://pedcbioportal.org/study/summary?id={study}#filterJson={encoded_filter}"

        with open("generated_url.txt", "w") as f:
            f.write(url + "\n")
        with open("generated_filter.json", "w") as f:
            json.dump(filter_json, f, indent=2)
        print("\nSaved URL to 'generated_url.txt' and filter JSON to 'generated_filter.json'")

    except (Exception, psycopg2.DatabaseError) as error:
        print(error)
        conn.close()
    conn.close()


def main():
    parser = argparse.ArgumentParser(description="Generate PedcBioPortal URL with embedded filters")
    parser.add_argument("-db", "--db-ini", action="store", dest="db_ini", help="Database config file - formatting like aws or sbg creds")
    parser.add_argument("-p", "--profile", action="store", dest="profile", help="ini profile name", default="postgresql")
    parser.add_argument("-s", "--study", action="store", dest="study", help="Cancer study ID to compare on server", default="pbta_all")
    
    args = parser.parse_args()
    run_py(args)


if __name__ == "__main__":
    main()
