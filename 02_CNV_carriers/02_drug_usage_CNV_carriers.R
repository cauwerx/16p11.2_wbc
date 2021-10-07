# Extract information on drug usage for 16p11.2 CNV carriers

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
### STEP 2: Extract medication columns ##########

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


#################################################
### STEP 3: Annotate medication #################

###### Load data ################################

# Description of UKBB medication usage field 20003 (code 4); Coding 4 is available on teh UKBB portal
ukb_med <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/03_covariates/data/raw/UKBB_20003_coding4_medication.tsv", col.names = c("UKBB_code", "UKBB_description")))

# Annotation of UKBB drugs with ATC 3rd level code; Available from Appleby et al., 2019 (https://www.researchsquare.com/article/rs-9729/v1)  
atc_ukb <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/03_covariates/data/raw/UKBB_20003_ATC_3rd_annotation_Appleby_2019.csv", select = c(1,3), col.names = c("UKBB_code", "ATC_code")))

# Description of ATC 3rd level; Generated from information on https://go.drugbank.com/atc
atc_med <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/03_covariates/data/raw/ATC_3rd_level_210929.csv", col.names = c("ATC_description", "ATC_code")))
atc_ukb <- left_join(atc_ukb, atc_med, by = "ATC_code")

###### Annotate ################################

MED <- left_join(MED, ukb_med, by = "UKBB_code") 
MED <- left_join(MED, atc_ukb, by = "UKBB_code") 


#################################################
### STEP 4: Medication use in CNV carriers ######

# Select medication used by CNV carriers
med_CNV <- MED[MED$eid %in% c(dup_carriers, del_carriers), ]
med_DUP <- MED[MED$eid %in% dup_carriers, ]
med_DEL <- MED[MED$eid %in% del_carriers, ]

# Add ATC code to the ATC description
med_CNV[which(!is.na(med_CNV$ATC_code)), "ATC_description"] <- paste0(med_CNV[which(!is.na(med_CNV$ATC_code)), "ATC_description"], " (", med_CNV[which(!is.na(med_CNV$ATC_code)), "ATC_code"], ")")
med_DUP[which(!is.na(med_DUP$ATC_code)), "ATC_description"] <- paste0(med_DUP[which(!is.na(med_DUP$ATC_code)), "ATC_description"], " (", med_DUP[which(!is.na(med_DUP$ATC_code)), "ATC_code"], ")")
med_DEL[which(!is.na(med_DEL$ATC_code)), "ATC_description"] <- paste0(med_DEL[which(!is.na(med_DEL$ATC_code)), "ATC_description"], " (", med_DEL[which(!is.na(med_DEL$ATC_code)), "ATC_code"], ")")

# Replace missing ATC descriptions (NA) by "No ATC"
med_CNV[which(is.na(med_CNV$ATC_code)), "ATC_description"] <- "No ATC"
med_DUP[which(is.na(med_DUP$ATC_code)), "ATC_description"] <- "No ATC"
med_DEL[which(is.na(med_DEL$ATC_code)), "ATC_description"] <- "No ATC"

# Select drug categories used by at least n CNV/DUP/DEL carriers
n <- 1 
cat_CNV <- names(sort(table(med_CNV$ATC_description)[table(med_CNV$ATC_description) >= n]))
print(paste0("Number of categories with at least 1 CNV carrier: ", length(cat_CNV)))
cat_DUP <- names(sort(table(med_DUP$ATC_description)[table(med_DUP$ATC_description) >= n]))
print(paste0("Number of categories with at least 1 DUP carrier: ", length(cat_DUP)))
cat_DEL <- names(sort(table(med_DEL$ATC_description)[table(med_DEL$ATC_description) >= n]))
print(paste0("Number of categories with at least 1 DEL carrier: ", length(cat_DEL)))


#################################################
### Save ########################################

# Drug categories used by CNV carriers
fwrite(data.frame(ATC_description = cat_CNV[cat_CNV != "No ATC"]), "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/data/final/drug_usage/categories/drug_categories_1_CNV_carriers.txt", row.names = F, col.names = F, quote = F, sep = "\t")
fwrite(data.frame(ATC_description = cat_DUP[cat_DUP != "No ATC"]), "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/data/final/drug_usage/categories/drug_categories_1_DUP_carriers.txt", row.names = F, col.names = F, quote = F, sep = "\t")
fwrite(data.frame(ATC_description = cat_DEL[cat_DEL != "No ATC"]), "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/data/final/drug_usage/categories/drug_categories_1_DEL_carriers.txt", row.names = F, col.names = F, quote = F, sep = "\t")

# Exact medication usage records of CNV carriers
fwrite(med_CNV, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/data/final/drug_usage/drug_usage_CNV_carriers.txt", row.names = F, col.names = T, quote = F, sep = "\t")
fwrite(med_DUP, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/data/final/drug_usage/drug_usage_DUP_carriers.txt", row.names = F, col.names = T, quote = F, sep = "\t")
fwrite(med_DEL, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/data/final/drug_usage/drug_usage_DEL_carriers.txt", row.names = F, col.names = T, quote = F, sep = "\t")
