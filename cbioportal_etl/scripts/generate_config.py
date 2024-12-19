"""
Generate a JSON config from a TSV based on the provided study.

Command:
python3 generate_config.py --tsv all_studies_config_values.tsv --study oligo_nation
"""
import csv
import json
import argparse
import ast

def parse_value(value, should_split=False):
    value = value.strip()
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


def update_nested_dict(result, key, sub_key, value):
    should_split = ((sub_key and sub_key.lower() in {"mutation", "mafs.kf"}) or (key and key.lower() in {"dl_file_type_list"}))
    current = result.setdefault(key, {})
    if sub_key:
        keys = sub_key.split('.')
        for part in keys[:-1]:
            current = current.setdefault(part, {})
        current[keys[-1]] = parse_value(value, should_split)
    else:
        result[key] = parse_value(value, should_split)


def generate_json(args):
    result = {}
    filtered_rows = []
    with open(args.config_tsv, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if (row['Study'] != args.study and row['Study'] != "static") or not row['Value'].strip():
                continue
            filtered_rows.append(row)
            update_nested_dict(result, row['Key'], row['Sub-Key'], row['Value'])
    
    filtered_tsv_file = f"{args.study}_config_values.tsv"
    with open(filtered_tsv_file, 'w', newline='') as tsv_file:
        writer = csv.DictWriter(tsv_file, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(filtered_rows)

    output_file = f"{args.study}_config.json"
    with open(output_file, 'w') as json_file:
        json.dump(result, json_file, indent=4)
    print(f"JSON file generated: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Generate JSON from TSV based on study")
    parser.add_argument("-v", "--config-tsv", action="store", dest="config_tsv", help="Path to the input TSV file")
    parser.add_argument("-s", "--study", action="store", dest="study", help="Cancer study ID to compare on server")
    
    args = parser.parse_args()
    generate_json(args)


if __name__ == "__main__":
    main()
