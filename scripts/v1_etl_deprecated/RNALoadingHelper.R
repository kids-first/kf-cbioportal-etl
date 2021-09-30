##############################################
#Purpose: Code to create RNA-Seq raw data per disease
#Date: 10/5/2018
#Author: Pichai Raman
##############################################


#Read in mapping file
#ExpressionFile 

init_res <- function(fpkmFileLoc=NULL, mapFileLoc=NULL, genesFileLoc=NULL)
{
	print(paste("Reading inputs for expression data"))
	res = read.table(file=fpkmFileLoc, header=TRUE, sep="\t", check.names = FALSE)
	cBioGenes <- levels(read.delim(genesFileLoc)[,2]);

	#Format RNA-Seq Rows
	res[,"median"] <- apply(res[3:ncol(res)], FUN=median, MARGIN=1)
	res <- res[order(-res[,"median"]),]
	res <- res[!duplicated(as.character(res[,2])),]
	rownames(res) <- res[,2];
	res <- res[-1:-2];
	res <- res[intersect(rownames(res), cBioGenes),]
	res <- res[-ncol(res)]

	#load mapping file and rename columns
	mappingFile <- read.delim(mapFileLoc, header=F, stringsAsFactors=F)
	rownames(mappingFile) <- mappingFile[,2];
	res <- res[, intersect(colnames(res), rownames(mappingFile))]
	colnames(res) <- mappingFile[colnames(res),1]
	return(res)

}

myZ <- function(x)
{
	x <- (x-mean(x))/sd(x);
	return(x);
}

createRNASeqExpMatrix <- function(loc=NULL, myDisease=NULL, dataSheetsDir=NULL, res=NULL)
{
	print(paste("Reading datasheet info for", myDisease, sep=" "))
	# z score out check
	writeOut <- 0
	#Filter and write out
	tmpData <- read.delim(paste(dataSheetsDir, myDisease, "/data_clinical_sample.txt", sep=""), skip=4);
	tmpSamps <- as.character(tmpData[,"SAMPLE_ID"]);
	tmpMat <- data.frame(res[,intersect(tmpSamps, colnames(res))]);
	if(ncol(tmpMat)==1) 
	{
		colnames(tmpMat)<- intersect(tmpSamps, colnames(res));
	  	rownames(tmpMat) <- rownames(res);
	}
	if(length(intersect(tmpSamps, colnames(res)))>0)
	{
		dataFPKM <- data.frame(rownames(tmpMat), tmpMat);

		#Write out expression matrix
		colnames(dataFPKM)[1] <- c("Hugo_Symbol")
		colnames(dataFPKM) <- gsub("^X", "", colnames(dataFPKM))
		colnames(dataFPKM) <- gsub("\\.", "-", colnames(dataFPKM))
		write.table(dataFPKM, paste(loc, "/data_rna_seq_v2_mrna.txt", sep=""), sep="\t", row.names=F, quote=F)

		#Write out Z-score

		if (length(intersect(tmpSamps, colnames(res))) > 1){
			tmpMat <- data.frame(log2(tmpMat+1));
			tmpMat <- data.frame(t(apply(tmpMat, FUN=myZ, MARGIN=1)))
			dataZ <- data.frame(rownames(tmpMat), tmpMat);
			colnames(dataZ) <- colnames(dataFPKM)
			write.table(dataZ, paste(loc, "/data_rna_seq_v2_mrna_median_Zscores.txt", sep=""), sep="\t", row.names=F, quote=F)
			writeOut <- 1
			print(paste("Finished writing out Expression Files for", myDisease, sep=" "));
		}
		else{
			print(paste("One or fewer samples for", myDisease, "detected, skipping Z score output!", sep=" "))
		}
	}
	return(writeOut)
}

#Create Meta File for RNA exp files
createExpressionMeta <- function(myLoc=NULL, cancStudyID=NULL, myDesc=NULL)
{
	#Create File
	myFile <- paste(myLoc, "/meta_rna_seq_v2_mrna.txt", sep="");
	file.create(myFile);

	#Add Cancer Study
	cancStudyIDLine <- paste("cancer_study_identifier: ", cancStudyID, sep="")
	write(cancStudyIDLine, myFile, append=TRUE)

	#Add stable id
	stableIDTypeLine <- "stable_id: rna_seq_mrna"
	write(stableIDTypeLine, myFile, append=TRUE)

	#Add Profile Name
	profileNameLine <- "profile_name: mRNA expression"
	write(profileNameLine, myFile, append=TRUE)

	#Add Profile Description
	profileDescLine <- myDesc
	write(profileDescLine, myFile, append=TRUE)

	#Add Genetic Alteration Type
	genAltTypeLine <- "genetic_alteration_type: MRNA_EXPRESSION";
	write(genAltTypeLine, myFile, append=TRUE)

	#Add Data Type
	dataTypeLine <- "datatype: CONTINUOUS"
	write(dataTypeLine, myFile, append=TRUE)

	#Add Show Profile in analysis Tab
	showProfileInTabLine <- "show_profile_in_analysis_tab: false"
	write(showProfileInTabLine, myFile, append=TRUE)

	#Add Data_filename
	dataFileLine <- "data_filename: data_rna_seq_v2_mrna.txt"
	write(dataFileLine, myFile, append=TRUE)

	print(paste("Done Creation RNA Expression Meta for", cancStudyID, sep=" "));
}

#Create Meta File for RNA exp files
createExpressionZMeta <- function(myLoc=NULL, cancStudyID=NULL, myDesc=NULL)
{
	#Create File
	myFile <- paste(myLoc, "/meta_rna_seq_v2_mrna_median_Zscores.txt", sep="");
	file.create(myFile);

	#Add Cancer Study
	cancStudyIDLine <- paste("cancer_study_identifier: ", cancStudyID, sep="")
	write(cancStudyIDLine, myFile, append=TRUE)

	#Add stable id
	stableIDTypeLine <- "stable_id: rna_seq_mrna_median_Zscores"
	write(stableIDTypeLine, myFile, append=TRUE)

	#Add Profile Name
	profileNameLine <- "profile_name: mRNA expression z-scores"
	write(profileNameLine, myFile, append=TRUE)

	#Add Profile Description
	profileDescLine <- myDesc;
	write(profileDescLine, myFile, append=TRUE)

	#Add Genetic Alteration Type
	genAltTypeLine <- "genetic_alteration_type: MRNA_EXPRESSION";
	write(genAltTypeLine, myFile, append=TRUE)

	#Add Data Type
	dataTypeLine <- "datatype: Z-SCORE"
	write(dataTypeLine, myFile, append=TRUE)

	#Add Show Profile in analysis Tab
	showProfileInTabLine <- "show_profile_in_analysis_tab: true"
	write(showProfileInTabLine, myFile, append=TRUE)

	#Add Data_filename
	dataFileLine <- "data_filename: data_rna_seq_v2_mrna_median_Zscores.txt"
	write(dataFileLine, myFile, append=TRUE)

	print(paste("Done Creation Expression Z Meta", cancStudyID, sep=" "));
}






