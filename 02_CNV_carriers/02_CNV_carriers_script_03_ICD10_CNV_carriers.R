# Extract information on ICD10 diagnoses for 16p11.2 CNV carriers

#################################################
### Libraries ###################################
library(data.table)
library(tidyr)
library(dplyr)


#################################################
### STEP 1: Load CNV carriers ###################

# Load retained sample eids; This file was genrated by 02_CNV_carriers/01_CNV_carriers
selected_eid <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/data/final/16p11.2_CNV_profile_HC.txt", header = T, select = c(1,2,5), col.names = c("eid", "sex", "CNV")))

# Select high confidence 16p11.2 duplication and deletion carriers
dup_carriers <- selected_eid[which(selected_eid$CNV == 1), "eid"]
del_carriers <- selected_eid[which(selected_eid$CNV == -1), "eid"]


#################################################
### STEP 2: Extract ICD10 columns ###############

# Identify columns with ICD10 diagnoses information (field 41270); "ukb44073.csv", a phenotype file available from the UKBB portal
header <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb44073.csv", header = T, nrow = 0))
col_ICD <- grep("^41270-", names(header))

# Extract ICD10 diagnoses columns
ICD <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb44073.csv", header = T, select = c(1, col_ICD), colClass = "character"))
ICD[ICD==""] <- NA
ICD$eid <- as.numeric(ICD$eid)

# Retain only CNV carriers rows
ICD_CNV <- ICD[ICD$eid %in% c(del_carriers, dup_carriers), ]

# Remodel
ICD_CNV <- na.omit(gather(ICD_CNV, Instance, Code, 2:ncol(ICD_CNV)))
ICD_CNV <- unique(ICD_CNV[, c(1,3)])


#################################################
### Save ########################################

fwrite(ICD_CNV, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/data/final/ICD10_diagnoses/ICD10_diagnoses_CNV_carriers.txt", col.names = F, row.names = F, quote = F, sep = "\t")


