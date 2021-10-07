# Extract drug usage information for all individuals to be used as a covariate

#################################################
### Libraries ###################################
library(data.table)
library(tidyr)
library(dplyr)


#################################################
### STEP 1: Load data ###########################

#### Load all samples ###########################

# This file contains all retained eids and was generate by 01_samples/01_sample_filtering.R
eids <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/01_samples/data/final/samples_white_british_All.txt", header = F, col.names = "eid"))


#### Load CNV ATC categories & Merge ############

# This file contains all ATC categories used by at least 1 CNV carrier and was generate by 02_CNV_carriers/02_drug_usageCNV_carriers.R
cat <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/data/final/drug_usage/categories/drug_categories_1_CNV_carriers.txt", sep = "\t", header = F, col.names = "cat"))


#### Load UKBB-ATC files & merge ################

# Annotation of UKBB drugs with ATC 3rd level code; Available from Appleby et al., 2019 (https://www.researchsquare.com/article/rs-9729/v1)   
atc_ukb <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/03_covariates/data/raw/UKBB_20003_ATC_3rd_annotation_Appleby_2019.csv", select = c(1,3), col.names = c("UKBB_code", "ATC_code")))

# Description of ATC 3rd level; Generated from information on https://go.drugbank.com/atc
atc_desc <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/03_covariates/data/raw/ATC_3rd_level_210929.csv", col.names = c("ATC_description", "ATC_code")))

# Merge
atc_ukb <- left_join(atc_ukb, atc_desc, by = "ATC_code")


#################################################
### STEP 2: Extract individual drug usage ######

# Identify columns with medication usage information (field 20003); "ukb44073.csv", a phenotype file available from the UKBB portal 
header <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb44073.csv", header = T, nrow = 0))
col_MED <- grep("^20003-", names(header))

# Extract these columns from the phenotype file and remodel the data; "ukb44073.csv", a phenotype file available from the UKBB portal  
MED <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb44073.csv", header = T, select = c(1, col_MED), colClass = "character"))
MED[MED==""] <- NA
MED <- na.omit(gather(MED, Instance, Code, 2:ncol(MED)))
MED <- unique(MED[, c(1,3)]) #  1'325'717 entries
MED <- MED %>% mutate_if(is.character, as.numeric)
colnames(MED) <- c("eid", "UKBB_code")

# Add ATC annotation
MED <- left_join(MED, atc_ukb[,c(1,2)], by = "UKBB_code")


#################################################
### STEP 3: Add ATC categories as a covariates ##

# The covariate file is generated from the eids file
cov <- eids 

# Loop over the different ATC categories used by CNV carriers
for (c in cat$cat) {
	
	# Define the ATC category
	print(paste0("Adding ", c))
	cat_num <- gsub(".*\\(", "", gsub(")", "",c))

	# Add an empty column for that category
	cov[, cat_num] <- 0

	# Identify individuals reporting drug usage from that category
	users <- unique(MED[which(MED$ATC_code == cat_num), "eid"])
	print(paste0("Number of drug users (full UK Biobank) in category", cat_num, ": ", length(users)))

	# Indicate drug usage in the covariate file
	cov[cov$eid %in% users, cat_num] <- 1
}


#################################################
### Save ########################################

fwrite(cov, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/03_covariates/data/final/drug_covariates.txt", row.names = F, col.names = T, quote = F, sep = "\t")
