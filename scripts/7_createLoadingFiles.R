#!/usr/bin/env Rscript
##############################################
#Purpose: Code to create Pedcbioportal loading files
#Date: 10/5/2018
#Author: Pichai Raman, Modified by Miguel Brown
##############################################

##############################################
#Rough outline
#1. We will start with the clinical file and update it with meta file
#2. We will find corresponding MAF files and create MAF meta file
#3. We will then find corresponding RNA-Seq files and update
##############################################

#Set Mutation, CNA, and Datasheets Directory
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Please supply a json config file\n", call.=FALSE)
}
library("rjson")
config_data  <- fromJSON(file=args[1])
data_dir = config_data$data_dir
mutDir <- paste(data_dir, "/merged_mafs/", sep = "")
cnaDir <- paste(data_dir, "/merged_cnvs/", sep = "")
hasRNA = config_data$rna_flag
if (hasRNA){
	fusion_dir <- paste(data_dir, config_data$fusion_dir, sep = "")
}
hasCNA = config_data$cna_flag
dataSheetsDir <- paste(data_dir, "/datasheets/", sep = "")
diseaseMappingFile <- paste(data_dir, "/", config_data$dx_tbl_fn, sep = "");


#Source RNA-Seq Functions
script_dir = config_data$script_dir
# source("RNALoadingHelper.R")
source(paste(script_dir,"/CNALoadingHelper.R", sep = ""))

#Create Meta File for MAF files
createMutationsMeta <- function(myLoc=NULL, cancStudyID=NULL, myDesc=NULL)
{
	#Create File
	myFile <- paste(myLoc, "/meta_mutations_extended.txt", sep="");
	file.create(myFile);

	#Add Cancer Study
	cancStudyIDLine <- paste("cancer_study_identifier: ", cancStudyID, sep="")
	write(cancStudyIDLine, myFile, append=TRUE)

	#Add stable id
	stableIDTypeLine <- "stable_id: mutations"
	write(stableIDTypeLine, myFile, append=TRUE)

	#Add Profile Name
	profileNameLine <- "profile_name: Mutations"
	write(profileNameLine, myFile, append=TRUE)

	#Add Profile Description
	profileDescLine <- paste("profile_description: ", myDesc, sep="")
	write(profileDescLine, myFile, append=TRUE)

	#Add Genetic Alteration Type
	genAltTypeLine <- "genetic_alteration_type: MUTATION_EXTENDED";
	write(genAltTypeLine, myFile, append=TRUE)

	#Add Data Type
	dataTypeLine <- "datatype: MAF"
	write(dataTypeLine, myFile, append=TRUE)

	#Add Show Profile in analysis Tab
	showProfileInTabLine <- "show_profile_in_analysis_tab: true"
	write(showProfileInTabLine, myFile, append=TRUE)

	#Add Data_filename
	dataFileLine <- "data_filename: data_mutations_extended.txt"
	write(dataFileLine, myFile, append=TRUE)

	print(paste("Done Creation Mutations Meta for", cancStudyID, sep = " "));
}

#Create Meta File for MAF files
createMetaStudy <- function(myLoc=NULL, myTypeOfCancer=NULL, cancStudyID=NULL, myName=NULL, myDesc=NULL, myShortName=NULL, myGroups=NULL)
{
	#Create File
	myFile <- paste(myLoc, "/meta_study.txt", sep="");
	file.create(myFile);

	#Add Type of Cancer
	typeOfCancerLine <- paste("type_of_cancer: ", myTypeOfCancer, sep="")
	write(typeOfCancerLine, myFile, append=TRUE)

	#Add Cancer Study
	cancStudyIDLine <- paste("cancer_study_identifier: ", cancStudyID, sep="")
	write(cancStudyIDLine, myFile, append=TRUE)

	#Add Profile Name
	nameLine <- paste("name: ", myName, sep="")
	write(nameLine, myFile, append=TRUE)

	#Add Profile Description
	descLine <- paste("description: ", myDesc, sep="")
	write(descLine, myFile, append=TRUE)

	#Add Data_filename
	shortNameLine <- paste("short_name: ", myShortName, sep="")
	write(shortNameLine, myFile, append=TRUE)

	#Add group
	groupLine <- paste("groups: ", myGroups, sep="")
	write(groupLine, myFile, append=TRUE)

	print(paste("Done Creation Meta Study for", cancStudyID, sep = " "));
}

#Create Meta Clinical Sample
createMetaSample <- function(myLoc=NULL, cancStudyID=NULL)
{
	#Create File
	myFile <- paste(myLoc, "/meta_clinical_sample.txt", sep="");
	file.create(myFile);

	#Add Cancer Study
	cancStudyIDLine <- paste("cancer_study_identifier: ", cancStudyID, sep="")
	write(cancStudyIDLine, myFile, append=TRUE)

	#Add Genetic Alteration type
	genAltLine <- "genetic_alteration_type: CLINICAL"
	write(genAltLine, myFile, append=TRUE)

	#Add Profile Description
	datTypeLine <- "datatype: SAMPLE_ATTRIBUTES"
	write(datTypeLine, myFile, append=TRUE)

	#Add Data_filename
	datFileNameLine <- "data_filename: data_clinical_sample.txt"
	write(datFileNameLine, myFile, append=TRUE)
	print(paste("Done Creation Meta Sample for", cancStudyID, sep=" "));
}

#Create Meta Clinical Patient
createMetaClinical <- function(myLoc=NULL, cancStudyID=NULL)
{
	#Create File
	myFile <- paste(myLoc, "/meta_clinical_patient.txt", sep="");
	file.create(myFile);

	#Add Cancer Study
	cancStudyIDLine <- paste("cancer_study_identifier: ", cancStudyID, sep="")
	write(cancStudyIDLine, myFile, append=TRUE)

	#Add Genetic Alteration type
	genAltLine <- "genetic_alteration_type: CLINICAL"
	write(genAltLine, myFile, append=TRUE)

	#Add Profile Description
	datTypeLine <- "datatype: PATIENT_ATTRIBUTES"
	write(datTypeLine, myFile, append=TRUE)

	#Add Data_filename
	datFileNameLine <- "data_filename: data_clinical_patient.txt"
	write(datFileNameLine, myFile, append=TRUE)
	print(paste("Done Creation Meta Patient for", cancStudyID, sep=" "));
}


createCaseLists <- function(myLoc=NULL, cancStudyID=NULL, cancStableID=NULL, caseListName=NULL, caseListDesc=NULL, caseListCat=NULL, myCaseListIDs=NULL)
{
	#Create File
	myFile <- paste(myLoc, sep="");
	file.create(myFile);

	#Add Cancer Study
	cancStudyIDLine <- paste("cancer_study_identifier: ", cancStudyID, sep="")
	write(cancStudyIDLine, myFile, append=TRUE)

	#Add Cancer Stable ID
	cancStableIDLine <- paste("stable_id: ", cancStableID, sep="")
	write(cancStableIDLine, myFile, append=TRUE)

	#Add Case List Name
	caseNameLine <- paste("case_list_name: ", caseListName, sep="")
	write(caseNameLine, myFile, append=TRUE)

	#Add Case List Description
	caseDescLine <- paste("case_list_description: ", caseListDesc, sep="")
	write(caseDescLine, myFile, append=TRUE)

	#Add Case List Category
	caseCatLine <- paste("case_list_category: ", caseListCat, sep="")
	write(caseCatLine, myFile, append=TRUE)

	#Add list ids
	myCaseListIDs <- paste(myCaseListIDs, collapse="\t")
	caseListIDsLine <- paste("case_list_ids: ", myCaseListIDs, sep="")
	write(caseListIDsLine, myFile, append=TRUE)

	print(paste("Done Creation of case list for", cancStudyID, sep=" "));

}



createStudyAll <- function(myDisease=NULL, 
	myName=NULL, 
	myDescMS=NULL, 
	myDescMM=config_data$mut_desc,
	myDescMRF="profile_description: Expression levels from RNA-Seq (FPKM)",
	myDescMRZ="profile_description: Expression levels from RNA-Seq (z-scores)",
	myDescMC="profile_description: Predicted copy number values from WGS (Continuous)",
	myDescMCD="profile_description: Predicted copy number values from WGS (Discrete)")
{
	myDiseaseStudyFolder <- paste("./processed/", config_data$cancStudyID, sep="")
	myCancStudyID <- config_data$cancStudyID;
	if (config_data$prepend_cancStudyID_flag){
		myDiseaseStudyFolder <- paste("./processed/", myDisease, "_", config_data$cancStudyID, sep="")
		myCancStudyID <- paste(myDisease, "_", config_data$cancStudyID, sep="");
	}
	caseListFolder <- paste(myDiseaseStudyFolder, "/case_lists", sep="")

	#Create Directory
	dir.create(myDiseaseStudyFolder, recursive = TRUE)
	dir.create(caseListFolder, recursive = TRUE);
	#Move clinical and patient files
	system(paste("cp ", dataSheetsDir, myDisease, "/data* ", myDiseaseStudyFolder, sep=""));

	#Move MAF File
	system(paste("cp ", mutDir, myDisease, ".strelka.vep.filtered.maf ", myDiseaseStudyFolder, "/data_mutations_extended.txt", sep=""));

	if(hasCNA){
		#Move CNA File
		system(paste("cp ", cnaDir, myDisease, ".predicted_cnv.txt ", myDiseaseStudyFolder, "/data_linear_CNA.txt", sep=""));
		system(paste("cp ", cnaDir, myDisease, ".discrete_cnvs.txt ", myDiseaseStudyFolder, "/data_CNA.txt", sep=""));

		#Create Discrete CNA File
		# tmpCNV <- read.delim(paste(cnaDir, myDisease, ".predicted_cnv.txt",sep=""))
		# tmpCNVGenes <- tmpCNV[1];
		# tmpCNVVals <- as.matrix(tmpCNV[2:ncol(tmpCNV)]);
		# tmpCNVVals <- tmpCNVVals-2;
		# tmpCNVVals[tmpCNVVals>6]<-2
		# tmpCNVVals[tmpCNVVals>2.1]<-1
		# tmpCNV <- data.frame(tmpCNVGenes, tmpCNVVals)
		# colnames(tmpCNV) <- gsub("^X", "", colnames(tmpCNV))
		# colnames(tmpCNV) <- gsub("\\.", "-", colnames(tmpCNV))
		# write.table(tmpCNV, paste(myDiseaseStudyFolder, "/data_CNA.txt", sep=""), row.names=F, sep="\t", quote=F)
		#Create CNV meta files
		createCNVMeta(myDiseaseStudyFolder, cancStudyID=myCancStudyID, myDescMC);
		createCNVDiscreteMeta(myDiseaseStudyFolder, cancStudyID=myCancStudyID, myDescMCD);
		#Case List CNA
		tmpData <- read.delim(paste(myDiseaseStudyFolder, "/data_linear_CNA.txt", sep=""));
		tmpSamps <- as.character(colnames(tmpData));
		tmpSamps <- tmpSamps[2:length(tmpSamps)]
		tmpSamps <- gsub("^X", "", tmpSamps)
		tmpSamps <- gsub("\\.", "-", tmpSamps)
		createCaseLists(paste(myDiseaseStudyFolder, "/case_lists/cases_cna.txt", sep=""),
			cancStudyID=myCancStudyID,
			cancStableID=paste(myCancStudyID, "_cna", sep=""),
			caseListName="Tumor Samples with CNA data",
			caseListDesc=paste("All tumors with CNA data (", length(tmpSamps), " samples)", sep=""),
			caseListCat="all_cases_with_cna_data",
			myCaseListIDs=tmpSamps)
		}

	#Create Study met file
	createMetaStudy(myDiseaseStudyFolder, 
		myTypeOfCancer=myDisease, 
		cancStudyID=myCancStudyID, 
		myName=myName, 
		myDesc=myDescMS, 
		myShortName=myCancStudyID,
		myGroups=config_data$group
		);

	#Create mutation meta file
	createMutationsMeta(myDiseaseStudyFolder, 
		cancStudyID=myCancStudyID, 
		myDesc=myDescMM);

	if(hasRNA)
	{
		#Create Expression meta files
		z_check = createRNASeqExpMatrix(loc=myDiseaseStudyFolder, myDisease=myDisease, dataSheetsDir=dataSheetsDir, res=res);
		createExpressionMeta(myDiseaseStudyFolder, cancStudyID=myCancStudyID, myDescMRF);
		if (z_check){
			createExpressionZMeta(myDiseaseStudyFolder, cancStudyID=myCancStudyID, myDescMRZ);
		}
		#Copy fusions files
		# system(paste("cp ", fusion_dir, myDisease, "/data_fusions.txt ", myDiseaseStudyFolder, "/data_fusions.txt", sep=""));
		# system(paste("cp ", fusion_dir, myDisease, "/meta_FUSION.txt ", myDiseaseStudyFolder, "/meta_FUSION.txt", sep=""));
	}


	#Create meta sample & Clinical
	createMetaSample(myDiseaseStudyFolder, cancStudyID=myCancStudyID);
	createMetaClinical(myDiseaseStudyFolder, cancStudyID=myCancStudyID);

	#Create case lists
	#Case List All
	tmpData <- read.delim(paste(myDiseaseStudyFolder, "/data_clinical_sample.txt", sep=""), skip=4);
	tmpSamps <- as.character(unique(tmpData[,"SAMPLE_ID"]));
	createCaseLists(paste(myDiseaseStudyFolder, "/case_lists/cases_all.txt", sep=""), 
		cancStudyID=myCancStudyID,
		cancStableID=paste(myCancStudyID, "_all", sep=""), 
		caseListName="All Tumors", 
		caseListDesc=paste("All tumor samples (", length(tmpSamps), " samples)", sep=""),
		caseListCat="all_cases_in_study",
		myCaseListIDs=tmpSamps)
	print("Done Creation cases all");

	#Case List Mutations
	print("Creating mutations case list");
	tmpData <- read.delim(paste(mutDir, myDisease, ".strelka.vep.filtered.maf", sep=""), skip=1, row.names=NULL);
	tmpSamps <- as.character(unique(tmpData[,"Tumor_Sample_Barcode"]));
	print("Read in mutations");
	#For now, fusions are lumped into mutation data, need to concat together
	if (hasRNA){
		print("RNA flag given, reading fusions as mutations");
		tmpData1 <- read.delim(paste(fusion_dir, myDisease, "/data_fusions.txt", sep=""), row.names=NULL);
		tmpSamps1 <- as.character(unique(tmpData[,"Tumor_Sample_Barcode"]));
		tmpSamps <- unique(rbind(tmpSamps, tmpSamps1))
		print("Done reading fusions");
	}
	createCaseLists(paste(myDiseaseStudyFolder, "/case_lists/cases_sequenced.txt", sep=""), 
		cancStudyID=myCancStudyID,
		cancStableID=paste(myCancStudyID, "_sequenced", sep=""), 
		caseListName="All Sequenced Tumors", 
		caseListDesc=paste("All sequenced tumor samples (", length(tmpSamps), " samples)", sep=""),
		caseListCat="all_cases_with_mutation_data",
		myCaseListIDs=tmpSamps)
	print("Created sequenced case list");

	if(hasRNA)
	{
	#Case List RNA
	print("Creating RNA case list");
	tmpData <- read.delim(paste(myDiseaseStudyFolder, "/data_rna_seq_v2_mrna.txt", sep=""));
	tmpSamps <- as.character(colnames(tmpData));
	tmpSamps <- tmpSamps[2:length(tmpSamps)]
	tmpSamps <- gsub("^X", "", tmpSamps)
	tmpSamps <- gsub("\\.", "-", tmpSamps)
	createCaseLists(paste(myDiseaseStudyFolder, "/case_lists/cases_RNA_Seq_v2_mRNA.txt", sep=""),
		cancStudyID=myCancStudyID,
		cancStableID=paste(myCancStudyID, "_rna_seq_v2_mrna", sep=""),
		caseListName="Tumor Samples with mRNA data (RNA Seq V2)",
		caseListDesc=paste("All samples with mRNA expression data (", length(tmpSamps), " samples)", sep=""),
		caseListCat="all_cases_with_mrna_rnaseq_data",
		myCaseListIDs=tmpSamps)
	}
	print("Created RNA case list");
}


#Create studies
diseaseMapping <- read.delim(diseaseMappingFile, header=T, stringsAsFactors=F)
colnames(diseaseMapping) <- c("Disease", "Abbr", "Disease2");
diseaseMapping <- diseaseMapping[!duplicated(diseaseMapping[,2]),]
rownames(diseaseMapping) <- diseaseMapping[,"Abbr"];
myFiles <- list.files(dataSheetsDir)

if(hasRNA)
{
  #Create Expression meta files
  source(paste(script_dir,"/RNALoadingHelper.R", sep = ""))
  res = init_res(config_data$fpkmFileLoc, config_data$mapFileLoc, config_data$genesFileLoc)
}
for(i in 1:length(myFiles))
{
	print(i);
	myAbbrev <- myFiles[i];
	myDisease <- as.character(diseaseMapping[myAbbrev,3]);
	createStudyAll(myDisease=myAbbrev, 
	myName=paste(myDisease, config_data$study_short_desc, sep=" "),
	myDescMS=paste(config_data$study_desc_p1, myDisease, config_data$study_desc_p2, sep=" ")
	)

}

