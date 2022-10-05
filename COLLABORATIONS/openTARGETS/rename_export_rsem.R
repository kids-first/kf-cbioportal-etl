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
    help = "mapping ID file with headers: BS_ID	Sample Type	Cbio ID"
  ),
  make_option(
    opt_str = "--type",
    type = "character",
    help = "study name, like 'openpbta'"
  )
)
#parse options
opts <- parse_args(OptionParser(option_list = option_list))
message("Reading in rds file ", opts$rna_rds)
rna = readRDS(opts$rna_rds)

# Example entry: brain	BS_6H1C1ME9	RNA	7316-2558_460366
message("Reading in sample naming file ", opts$map_id)
map_ids = read.csv(sep = "\t", header=TRUE, opts$map_id)
map_ids = map_ids[map_ids$Sample.Type == "RNA",]
# Subset only on sample present in mapping ID file
subset_rna = rna[,which(colnames(rna) %in% map_ids$BS_ID)]
rownames(subset_rna) = rownames(rna)
# Rename BS IDs to cBio IDs
message("Renaming samples and outputting tpm file")
setnames(subset_rna, old=as.character(map_ids$BS_ID), new=as.character(map_ids$Cbio.ID), skip_absent=TRUE)

write_tsv(data.frame("Hugo_Symbol"=rownames(subset_rna),subset_rna, check.names = FALSE),paste(opts$type, ".rsem_merged.txt", sep=""), escape="none")
# Get z score of log2 tpm with added pseudocount - round to 4 places as added precision not needed

rm(rna)
rm(opts)
rm(map_ids)

message("Calculating z scores of log2 tmp + 1")
#subset_zscore = round(t(scale(t(log2(subset_rna + 1)))), 4)
genes=rownames(subset_rna)

library("Rcpp")
cppFunction("
List compute_zscore(List num) {
int n = num.size();
float a =0;
float sum =0;
float st = 0;
float var = 0;
float mean =0;
for(int i = 0; i < n; ++i) {
  a=num[i];
  a = std::log2(a + 1);
  sum=sum+a;
  num[i]=a;
} 
mean=sum/n;
for(int i =0; i < n; ++i){
  a=num[i];
  var=var+(a-mean)*(a-mean);
}
st=std::sqrt(var/(n-1));
int round_off = 10000;
for(int i = 0; i < n; ++i) {
    a=num[i];
    a = std::ceil((a - mean)*round_off/st) ;
    num[i]=a/round_off;
}  
return num;
}")
zscore <- list()
test = t(subset_rna)
message(test[0])
stop()
#a=1
#message(Sys.time())
for (x in genes){
#message((a))
#subset_zscore=compute_zscore(subset_rna[x,]) # less than 3 mins for 500 genes C++ function
subset_zscore = round(scale(log2(subset_rna[x,] + 1)), 4) # less than 10 mins for 85 genes using scale function
#zscore <-append(zscore, subset_zscore)
#a=a+1
}
subset_zscore = t(subset_zscore)

stop("Stopping msg")
message("Outputting z scores")
write_tsv(data.frame("Hugo_Symbol"=rownames(subset_zscore),subset_zscore, check.names = FALSE),paste(opts$type, ".rsem_merged_zscore.txt", sep=""), escape="none")