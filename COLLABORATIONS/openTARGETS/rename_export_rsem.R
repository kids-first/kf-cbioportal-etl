if(!require(data.table)){install.packages('data.table')}
if(!require(readr)){install.packages('readr')}
rna = readRDS("../INPUTS/gene-expression-rsem-tpm-collapsed.rds")
# Read in table with header: Cbio_project	T_CL_BS_ID	Sample Type	Cbio_Tumor_Name	File_Type
# Example entry: brain	BS_6H1C1ME9	RNA	7316-2558_460366	rsem
# should have only RNA names in it
map_ids = read.csv(sep = "\t", header=TRUE, "../INPUTS/cbio_pbta_style_names.txt")
# For this load, subset only on pbta entries
pbta_rna = rna[,which(colnames(rna) %in% map_ids$T_CL_BS_ID)]
rownames(pbta_rna) = rownames(rna)
# Rename BS IDs to cBio IDs
setnames(pbta_rna, old=as.character(map_ids$T_CL_BS_ID), new=as.character(map_ids$Cbio_Tumor_Name))
write_tsv(data.frame("Hugo_Symbol"=rownames(pbta_rna),pbta_rna, check.names = FALSE),"data_rna_seq_v2_mrna.txt", quote_escape="none")
# Get z score of log2 tpm with added pseudocount - round to 4 places as added precision not needed
pbta_zscore = round(t(scale(t(log2(pbta_rna + 1)))), 4)

write_tsv(data.frame("Hugo_Symbol"=rownames(pbta_zscore),pbta_zscore, check.names = FALSE),"data_rna_seq_v2_mrna_median_Zscores.txt", quote_escape="none")