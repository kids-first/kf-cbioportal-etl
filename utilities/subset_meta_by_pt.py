import os
import sys


def subset_file(in_file: str, out_file: str, num_header_lines: int, check_list: list, check_index: int, add_index: int) -> list:
    """Subset files based on a list of values and a column to search.

    Output results to a new file and returns a list of values from another column in the matching rows.

    Args:
        in_file (str): Input file path.
        out_file (str): Output file path.
        num_header_lines (int): Number of header lines to copy to the output file.
        check_list (list): List of values to check against a specific column.
        check_index (int): Index of the column to check for values in check_list.
        add_index (int): Index of the column from which to extract values for the return list.

    Returns:
        list: A list of values from the specified add_index column in the matching rows.

    """
    add_list = []
    with open(in_file) as in_f, open(out_file, "w") as out_f:
        for _i in range(num_header_lines):
            head = next(in_f)
            out_f.write(head)
        for data in in_f:
            fields = data.split("\t")
            if fields[check_index] in check_list:
                add_list.extend(fields[add_index].split(";"))
                out_f.write(data)
    return add_list

cid_list = []
with open(sys.argv[1]) as cids:
    skip = next(cids)
    for line in cids:
        (cid, sid) = line.rstrip().split('\t')
        cid_list.append(cid)

datasheet_dir = sys.argv[2]
os.makedirs("subset/datasheets", exist_ok=True)
# data clinical patient and sample are always there
pt_list = subset_file(os.path.join(datasheet_dir, "data_clinical_patient.txt"), "subset/datasheets/data_clinical_patient.txt", 5, cid_list, 1, 0)
bs_list = subset_file(os.path.join(datasheet_dir, "data_clinical_sample.txt"), "subset/datasheets/data_clinical_sample.txt", 5, pt_list, 0, 3)
subset_file(sys.argv[3], "subset/cbio_file_name_id.txt", 1, bs_list, 1, 1)

# check datasets dir for all that are non-sample/patient and subset those too
for fname in os.listdir(datasheet_dir):
    if fname.startswith("data_") and fname.endswith(".txt") and fname not in ["data_clinical_patient.txt", "data_clinical_sample.txt"]:
        subset_file(os.path.join(datasheet_dir, fname), os.path.join("subset/datasheets", fname), 1, pt_list, 0, 0)
