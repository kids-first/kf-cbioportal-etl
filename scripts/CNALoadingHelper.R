##############################################
#Purpose: Code to create RNA-Seq raw data per disease
#Date: 9/5/2018
#Author: Pichai Raman
##############################################

#Create Meta File for MAF files
createCNVMeta <- function(myLoc=NULL, cancStudyID=NULL, myDesc=NULL)
{
	#Create File
	myFile <- paste(myLoc, "/meta_linear_CNA.txt", sep="");
	file.create(myFile);

	#Add Cancer Study
	cancStudyIDLine <- paste("cancer_study_identifier: ", cancStudyID, sep="")
	write(cancStudyIDLine, myFile, append=TRUE)

	#Add stable id
	stableIDTypeLine <- "stable_id: linear_CNA"
	write(stableIDTypeLine, myFile, append=TRUE)

	#Add Profile Name
	profileNameLine <- "profile_name: copy-number values"
	write(profileNameLine, myFile, append=TRUE)

	#Add Profile Description
	profileDescLine <- myDesc
	write(profileDescLine, myFile, append=TRUE)

	#Add Genetic Alteration Type
	genAltTypeLine <- "genetic_alteration_type: COPY_NUMBER_ALTERATION";
	write(genAltTypeLine, myFile, append=TRUE)

	#Add Data Type
	dataTypeLine <- "datatype: CONTINUOUS"
	write(dataTypeLine, myFile, append=TRUE)

	#Add Show Profile in analysis Tab
	showProfileInTabLine <- "show_profile_in_analysis_tab: false"
	write(showProfileInTabLine, myFile, append=TRUE)

	#Add Data_filename
	dataFileLine <- "data_filename: data_linear_CNA.txt"
	write(dataFileLine, myFile, append=TRUE)

	print("Done Creation CNV Meta");
}

#Create Meta File for MAF files
createCNVDiscreteMeta <- function(myLoc=NULL, cancStudyID=NULL, myDesc=NULL)
{
	#Create File
	myFile <- paste(myLoc, "/meta_CNA.txt", sep="");
	file.create(myFile);

	#Add Cancer Study
	cancStudyIDLine <- paste("cancer_study_identifier: ", cancStudyID, sep="")
	write(cancStudyIDLine, myFile, append=TRUE)

	#Add stable id
	stableIDTypeLine <- "stable_id: cna"
	write(stableIDTypeLine, myFile, append=TRUE)

	#Add Profile Name
	profileNameLine <- "profile_name: Binned copy-number values"
	write(profileNameLine, myFile, append=TRUE)

	#Add Profile Description
	profileDescLine <- myDesc
	write(profileDescLine, myFile, append=TRUE)

	#Add Genetic Alteration Type
	genAltTypeLine <- "genetic_alteration_type: COPY_NUMBER_ALTERATION";
	write(genAltTypeLine, myFile, append=TRUE)

	#Add Data Type
	dataTypeLine <- "datatype: DISCRETE"
	write(dataTypeLine, myFile, append=TRUE)

	#Add Show Profile in analysis Tab
	showProfileInTabLine <- "show_profile_in_analysis_tab: true"
	write(showProfileInTabLine, myFile, append=TRUE)

	#Add Data_filename
	dataFileLine <- "data_filename: data_CNA.txt"
	write(dataFileLine, myFile, append=TRUE)

	print("Done Creation CNV Meta Discrete");
}
