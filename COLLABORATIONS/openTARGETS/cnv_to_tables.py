import sys
import argparse
import gzip

parser = argparse.ArgumentParser(description='Script to convert openX cnv table to cBio format')
parser.add_argument('-m', '--mapping-file', action='store', dest = 'mapping_file', help='tsv file with header and bs_id, sample type, cbio ID mappings')
parser.add_argument('-c', '--copy-number', action='store', dest='cnv_tbl',
                    help='openX table')


def populate_id_map(map_fn):
    """
    Use the mapping file to create a BS ID to cBio ID mapping dict
    """
    map_file = open(map_fn)
    head = next(map_file)
    header = head.rstrip('\n').split('\t')
    b_idx = header.index('BS_ID')
    c_idx = header.index('Cbio ID')
    id_dict = {}
    for line in map_file:
        entry = line.rstrip('\n').split('\t')
        id_dict[entry[b_idx]] = entry[c_idx]
    map_file.close()
    return id_dict



def collate_data(cnv_fn):
    """
    Read in the cnv table, convert and hash vaules that are in the mapping dict
    """
    cnv_tbl = gzip.open(args.cnv_tbl)
    head = next(cnv_tbl).decode()
    header = head.rstrip('\n').split('\t')
    b_idx = header.index('biospecimen_id')
    s_idx = header.index('status')
    c_idx = header.index('copy_number')
    p_idx = header.index('ploidy')
    g_idx = header.index('gene_symbol')
    (c_dict, g_dict, p_dict) = ({}, {}, {})

    for line in cnv_tbl:
        data = line.decode().rstrip('\n').split('\t')
        if data[b_idx] in map_dict:
            samp_id = map_dict[data[b_idx]]
            gene = data[g_idx]
            ploidy = data[p_idx]
            cn = data[c_idx]
            gistic = qual_to_gistic[data[s_idx]]
            if samp_id not in p_dict:
                p_dict[samp_id] = ploidy
            if gene not in c_dict:
                c_dict[gene] = {}
                g_dict[gene] = {}
            c_dict[gene][samp_id] = cn
            g_dict[gene][samp_id] = gistic
    return (c_dict, g_dict, p_dict)


args = parser.parse_args()
sys.stderr.write("Getting mapping IDs\n")
map_dict = populate_id_map(args.mapping_file)

qual_to_gistic = {"deep deletion": "-2", "loss": "-1", "gain": "1", "amplification": "2"}
sys.stderr.write("Collating copy number data\n")
(cn_dict, gistic_dict, ploidy_dict) = collate_data(args.cnv_tbl)

sys.stderr.write("Printing results\n")
cnv_cn_out = open("data_linear_CNA.txt", "w")
cnv_cn_out.write("Hugo_Symbol")
cnv_gistic_out = open("data_CNA.txt", "w")
cnv_gistic_out.write("Hugo_Symbol")

sample_list = list(ploidy_dict.keys())
cnv_cn_out.write("\t" + "\t".join(sample_list) + "\n")
cnv_gistic_out.write("\t" + "\t".join(sample_list) + "\n")

# output results as gene x samples table
for gene in cn_dict:
    cnv_cn_out.write(gene)
    cnv_gistic_out.write(gene)
    for samp_id in sample_list:
        if samp_id in cn_dict[gene]:
            cnv_cn_out.write("\t" + cn_dict[gene][samp_id])
            cnv_gistic_out.write("\t" + gistic_dict[gene][samp_id])
        else:
            # if there is no CN values for that gene, assume neutral - 0 for gistic, equal to ploidy for raw CN
            cnv_cn_out.write("\t" + ploidy_dict[samp_id])
            cnv_gistic_out.write("\t0")
    cnv_cn_out.write("\n")
    cnv_gistic_out.write("\n")
cnv_cn_out.close()
cnv_gistic_out.close()