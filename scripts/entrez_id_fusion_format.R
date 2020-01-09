suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(org.Hs.eg.db))
hs <- org.Hs.eg.db


option_list <- list(
  make_option(c("-f", "--fusionfile"),type="character",
              help="Merged fusion calls from [STARfusion | Arriba] QC filtered for pedcbio"),
  make_option(c("-a", "--approvedHugoSymbols"), type="character",
              help="Hugo symbol approved gene name list"),
  make_option(c("-m","--manifest"),type="character",
               help="manifest with Tumor_Sample_Barcode\tCenter"),
  make_option(c("-o","--outputfile"),type="character",
              help="Standardized fusion calls from [STARfusion | Arriba] (.TSV)")
)

# Get command line options, if help option encountered print help and exit,
opt <- parse_args(OptionParser(option_list=option_list))
inputfile <- opt$fusionfile
approved <- opt$approvedHugoSymbols
manifest<-opt$manifest
outputfile <- opt$outputfile



total<-read_tsv(inputfile) %>% as.data.frame()
genes<-c(total$Gene1A,total$Gene1B,total$Gene2A,total$Gene2B)
manifest<-read_tsv(manifest) %>% as.data.frame() %>% unique()


my.symbols<-genes
entrez<-select(hs, keys=my.symbols, columns= c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
nomatch<-entrez[is.na(entrez$ENTREZID),]
#genes without entrez id
#nomatch

match<-read_tsv(approved) %>% as.data.frame()


my.symbols <- total$Gene1A
if (!all(is.na(my.symbols))){
  Gene1A_entrez<-select(hs, 
                      keys = my.symbols,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
  Gene1A_entrez[is.na(Gene1A_entrez$ENTREZID),"ENTREZID"]=""
} else {Gene1A_entrez=data.frame("ENTREZID"="","SYMBOL"="",stringsAsFactors = FALSE)}

my.symbols <- total$Gene2A
if (!all(is.na(my.symbols))){
  Gene2A_entrez<-select(hs,
                      keys = my.symbols,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
  Gene2A_entrez[is.na(Gene2A_entrez$ENTREZID),"ENTREZID"]=""
} else {Gene2A_entrez=data.frame("ENTREZID"="","SYMBOL"="",stringsAsFactors = FALSE)}

my.symbols <- total$Gene1B
if (!all(is.na(my.symbols))){
  Gene1B_entrez<-select(hs,
                      keys = my.symbols,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
  Gene1B_entrez[is.na(Gene1B_entrez$ENTREZID),"ENTREZID"]=""
} else {Gene1B_entrez=data.frame("ENTREZID"="","SYMBOL"="",stringsAsFactors = FALSE)}

my.symbols <- total$Gene2B
if (!all(is.na(my.symbols))){
  Gene2B_entrez<-select(hs,
                      keys = my.symbols,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")

  Gene2B_entrez[is.na(Gene2B_entrez$ENTREZID),"ENTREZID"]=""
}else {Gene2B_entrez=data.frame("ENTREZID"="","SYMBOL"="",stringsAsFactors = FALSE)}


# total$Gene1A_entrez<-Gene1A_entrez$ENTREZID
# total$Gene1B_entrez<-Gene1B_entrez$ENTREZID
# total$Gene2A_entrez<-Gene2A_entrez$ENTREZID
# total$Gene2B_entrez<-Gene2B_entrez$ENTREZID

total$Gene1A[is.na(total$Gene1A)]<-""
total$Gene2A[is.na(total$Gene2A)]<-""
total$Gene1B[is.na(total$Gene1B)]<-""
total$Gene2B[is.na(total$Gene2B)]<-""

total<-total %>% 
  left_join(unique(Gene1A_entrez), by=c("Gene1A"="SYMBOL")) %>% rename("ENTREZID"="Gene1A_entrez") %>%
  left_join(unique(Gene2A_entrez), by=c("Gene2A"="SYMBOL")) %>% rename("ENTREZID"="Gene2A_entrez") %>%
  left_join(unique(Gene1B_entrez), by=c("Gene1B"="SYMBOL")) %>% rename("ENTREZID"="Gene1B_entrez") %>%
  left_join(unique(Gene2B_entrez), by=c("Gene2B"="SYMBOL")) %>% rename("ENTREZID"="Gene2B_entrez")
  



#aggregate
fusion_format<-data.frame("Hugo_Symbol"=total$Gene1A,"Entrez_Gene_Id"=total$Gene1A_entrez,"Tumor_Sample_Barcode"=total$Sample,"Fusion"=total$FusionName,"DNA_support"=rep("no",nrow(total)),"RNA_support"=rep("yes",nrow(total)),"Method"=total$Caller,"Frame"=total$Fusion_Type)
fusion_format<-rbind(fusion_format,data.frame("Hugo_Symbol"=total$Gene2A,"Entrez_Gene_Id"=total$Gene2A_entrez,"Tumor_Sample_Barcode"=total$Sample,"Fusion"=total$FusionName,"DNA_support"=rep("no",nrow(total)),"RNA_support"=rep("yes",nrow(total)),"Method"=total$Caller,"Frame"=total$Fusion_Type))
fusion_format<-rbind(fusion_format,data.frame("Hugo_Symbol"=total$Gene1B,"Entrez_Gene_Id"=total$Gene1B_entrez,"Tumor_Sample_Barcode"=total$Sample,"Fusion"=total$FusionName,"DNA_support"=rep("no",nrow(total)),"RNA_support"=rep("yes",nrow(total)),"Method"=total$Caller,"Frame"=total$Fusion_Type))
fusion_format<-rbind(fusion_format,data.frame("Hugo_Symbol"=total$Gene2B,"Entrez_Gene_Id"=total$Gene2B_entrez,"Tumor_Sample_Barcode"=total$Sample,"Fusion"=total$FusionName,"DNA_support"=rep("no",nrow(total)),"RNA_support"=rep("yes",nrow(total)),"Method"=total$Caller,"Frame"=total$Fusion_Type))


fusion_format<-fusion_format %>% left_join(manifest,by=c("Tumor_Sample_Barcode"="Tumor_Sample_Barcode"))

fusion_format<-fusion_format[!is.na(fusion_format$Hugo_Symbol),]
fusion_format<-unique(fusion_format)

#filter other
if (nrow(fusion_format[grep("other",fusion_format$Frame),])>0){
fusion_format<-fusion_format[-grep("other",fusion_format$Frame),]
}
#filter intergenic
if (nrow(fusion_format[grep("/",fusion_format$Fusion),])>0){
fusion_format<-fusion_format[-grep("/", fusion_format$Fusion),]
}


#format for entrez== NA
fusion_format[is.na(fusion_format$Entrez_Gene_Id),"Entrez_Gene_Id"]<-""
fusion_format<-unique(fusion_format)

#remove blank Hugo Symbol rows caused by Gene2A/2B being empty for STARfusion and geneic arriba calls
fusion_format<-fusion_format[-which(fusion_format$Hugo_Symbol==""),]

#ensure order of cols correct before print
col_order = c("Hugo_Symbol", "Entrez_Gene_Id", "Center", "Tumor_Sample_Barcode", "Fusion", "DNA_support", "RNA_support", "Method", "Frame")
fusion_format = fusion_format[, col_order]
write.table(fusion_format,outputfile,quote = FALSE,row.names = FALSE,sep="\t")

# #get Center 
# #TGEN
# fusion_format[grep("C0",fusion_format$Tumor_Sample_Barcode),"Center"]<-"TGEN"
# #BGI
# #from Yuankun "all the fq source input are from BGI"
# manifest<-read.delim("~/Downloads/1562961569658-manifest.csv",sep=",",stringsAsFactors=F)
# manifest$Tumor_Sample_Barcode<-sub("_1.fq.gz","",manifest$name)
# manifest$Tumor_Sample_Barcode<-sub("_2.fq.gz","",manifest$Tumor_Sample_Barcode)
# manifest$Tumor_Sample_Barcode<-sub("_RNA","",manifest$Tumor_Sample_Barcode)
# manifest$Tumor_Sample_Barcode<-sub("7316_","7316-",manifest$Tumor_Sample_Barcode)
# fusion_format[which(fusion_format$Tumor_Sample_Barcode %in% manifest$Tumor_Sample_Barcode),"Center"]<-"BGI"
# #NANT
# fusion_format[-which(fusion_format$Center %in% c("TGEN","BGI")),"Center"]<-"NANT"


