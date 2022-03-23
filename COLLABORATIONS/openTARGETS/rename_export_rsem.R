if(!require(data.table)){install.packages('data.table')}
if(!require(readr)){install.packages('readr')}
if(!require(optparse)){install.packages('optparse')}
library("optparse")
library("data.table")

#process inputs
option_list <- list(
  make_option(
    opt_str = "--rna_rds",
    type = "character",
    help = "openX rsem rds expression object",
  ),
  make_option(
    opt_str = "--map_id",
    type = "character",
    help = "mapping ID file with headers: Cbio_project	BS_ID	Sample Type	Cbio.ID	File_Type"
  )
)
#parse options
opts <- parse_args(OptionParser(option_list = option_list))
message("Reading in rds file ", opts$rna_rds)
rna = readRDS(opts$rna_rds)
# Read in table with header: Cbio_project	BS_ID	Sample Type	Cbio.ID	File_Type
# Example entry: brain	BS_6H1C1ME9	RNA	7316-2558_460366	rsem
message("Reading in sample naming file ", opts$map_id)
map_ids = read.csv(sep = "\t", header=TRUE, opts$map_id)
map_ids = map_ids[map_ids$Sample.Type == "RNA",]
# Subset only on sample present in mapping ID file
subset_rna = rna[,which(colnames(rna) %in% map_ids$BS_ID)]
rownames(subset_rna) = rownames(rna)
# Rename BS IDs to cBio IDs
message("Renaming samples and outputting tpm file")
setnames(subset_rna, old=as.character(map_ids$BS_ID), new=as.character(map_ids$Cbio.ID))
write_tsv(data.frame("Hugo_Symbol"=rownames(subset_rna),subset_rna, check.names = FALSE),"mixed.rsem_merged.txt", escape="none")
# Get z score of log2 tpm with added pseudocount - round to 4 places as added precision not needed
message("Calculating z scores of log2 tmp + 1")
subset_zscore = round(t(scale(t(log2(subset_rna + 1)))), 4)
message("Outputting z scores")
write_tsv(data.frame("Hugo_Symbol"=rownames(subset_zscore),subset_zscore, check.names = FALSE),"mixed.rsem_merged_zscore.txt", escape="none")