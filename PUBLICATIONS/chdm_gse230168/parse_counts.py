import csv
import sys

from biomart import BiomartServer

server = BiomartServer("http://useast.ensembl.org/biomart")
print("Connected to server", file=sys.stderr)
fname = sys.argv[1]

# read in counts table
with open(fname, newline="\n") as counts_csv:
    creader = csv.reader(counts_csv, delimiter=",")
    header = next(creader)
    header[0] = "Hugo_Symbol"
    data = []
    gene_ids = []
    for row in creader:
        gene_ids.append(row[0].split(".")[0])
        data.append(row[1:])
print("read in csv", file=sys.stderr)
# get all id - sym mappings
attrs = ["ensembl_gene_id", "hgnc_symbol"]
mart = server.datasets["hsapiens_gene_ensembl"]
print("Got datasets from server", file=sys.stderr)
response = mart.search({"attributes": attrs})
print("Got ID-Sym pairs from server", file=sys.stderr)
check = response.raw.data.decode("ascii")
gene_id_sym_dict = {}

pairs_list = check.split("\n")
if pairs_list[-1] == "":
    del pairs_list[-1]
for x in pairs_list:
    try:
        (id, sym) = x.split("\t")
        gene_id_sym_dict[id] = sym
    except Exception as e:
        print(f"Failed on {x}", file=sys.stderr)

with open("counts_w_gene_sym.tsv", "w") as out:
    print("\t".join(header), file=out)
    for i in range(len(gene_ids)):
        if gene_ids[i] in gene_id_sym_dict and gene_id_sym_dict[gene_ids[i]] != "":
            sym = gene_id_sym_dict[gene_ids[i]]
        else:
            sym = gene_ids[i]
            print(f"Had to use ID:{sym}", file=sys.stderr)
        print(f"{sym}\t{'\t'.join(data[i])}", file=out)
