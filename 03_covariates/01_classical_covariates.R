# Extract classical covariates: sex, age, age^2, genotyping batch, PC1-40

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)


#################################################
### Load main data ##############################

# Sample QC; This corresponds to the "ukb_sqc_v2.txt" file available from the UKBB portal
sample_qc <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/geno/ukb_sqc_v2.txt", header = F, select = c(4,26:65), col.names = c("batch", paste0("PC", seq(1,40)))))

# Sample eid; This corresponds to a list of all sample eids (with sex information), retrieved from a ".fam" file available from the UKBB portal
sample_eid <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/plink/_001_ukb_cal_chr1_v2.fam", header = F, select = c(1,5), col.names = c("eid", "sex")))

# Merge
cov <- cbind(sample_eid, sample_qc)

# Selected retained sample eids; This file was genrated by 02_CNV_carriers/01_CNV_carriers
selected_eid <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/data/final/16p11.2_CNV_profile.txt", header = T, select = c(1), col.names = c("eid")))
cov <- cov[cov$eid %in% selected_eid$eid, ]


#################################################
### Add age and age^2 ###########################

# Identify column containing age information (field 21003) from "ukb28603.csv", a phenotype file available from the UKBB portal 
header <- fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb28603.csv", header = T, nrow = 0)
col_age <- grep("21003-", names(header))

# Extract age and determine age^2
age <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb28603.csv", header = T, select = c(1, col_age[1]), col.names = c("eid", "age")))
age$age2 <- age$age^2

# Add to main covariate file
cov <- right_join(age, cov, by = "eid")


#################################################
### Save ########################################

print(paste0("final dimensions: ", ncol(cov), " columns x ", nrow(cov), " rows"))
fwrite(cov, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/03_covariates/data/final/classical_covariates.txt", col.names = T, row.names = F, quote = F, sep = "\t")
