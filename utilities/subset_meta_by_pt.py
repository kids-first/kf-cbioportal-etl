"""Takes an input list if CIDs and external IDs, and subsets the datasheet files to only those CIDs and their associated samples. Outputs a new set of datasheet files in a subset/ directory.

Usage: python subset_meta_by_pt.py <cid_list> <datasheet_dir> <cbio_file_name_id>
"""

import os
import sys


def subset_file(in_file: str, out_file: str, num_header_lines: int, check_set: set, check_index: int, add_index: int) -> set:
    """Subset files based on a list of values and a column to search.

    Output results to a new file and returns a set of values from another column in the matching rows.

    Args:
        in_file (str): Input file path.
        out_file (str): Output file path.
        num_header_lines (int): Number of header lines to copy to the output file.
        check_set (list): List of values to check against a specific column.
        check_index (int): Index of the column to check for values in check_set.
        add_index (int): Index of the column from which to extract values for the return set.

    Returns:
        set: A set of values from the specified add_index column in the matching rows.

    """
    add_set: set[str] = set()
    with open(in_file) as in_f, open(out_file, "w") as out_f:
        for _i in range(num_header_lines):
            head = next(in_f)
            out_f.write(head)
        for data in in_f:
            fields = data.split("\t")
            if fields[check_index] in check_set:
                add_set.update(fields[add_index].split(";"))
                out_f.write(data)
    return add_set


def main():
    cid_set: set[str] = set()
    with open(sys.argv[1]) as cids:
        _skip = next(cids)
        for line in cids:
            (cid, _sid) = line.rstrip().split("\t")
            cid_set.add(cid)

    datasheet_dir = sys.argv[2]
    os.makedirs("subset/datasheets", exist_ok=True)
    # data clinical patient and sample are always there
    pt_set = subset_file(os.path.join(datasheet_dir, "data_clinical_patient.txt"), "subset/datasheets/data_clinical_patient.txt", 5, cid_set, 1, 0)
    bs_set = subset_file(os.path.join(datasheet_dir, "data_clinical_sample.txt"), "subset/datasheets/data_clinical_sample.txt", 5, pt_set, 0, 3)
    subset_file(sys.argv[3], "subset/cbio_file_name_id.txt", 1, bs_set, 1, 1)

    # check datasets dir for all that are non-sample/patient and subset those too
    for fname in os.listdir(datasheet_dir):
        if fname.startswith("data_") and fname.endswith(".txt") and fname not in ["data_clinical_patient.txt", "data_clinical_sample.txt"]:
            subset_file(os.path.join(datasheet_dir, fname), os.path.join("subset/datasheets", fname), 1, pt_set, 0, 0)

if __name__ == "__main__":
    main()
