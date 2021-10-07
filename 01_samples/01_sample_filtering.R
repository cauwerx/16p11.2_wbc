# Retrieve unrelated samples of white British ancestry

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### Load Main Data ##############################

# Sample QC; This corresponds to the "ukb_sqc_v2.txt" file available from the UKBB portal
sample_qc <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/geno/ukb_sqc_v2.txt", header = F, select = c(3:68)))
colnames(sample_qc) <- c("array", "batch", "plate", "well", "cluster_CR", "dQC", "dna_concentration", "submitted_gender", "inferred_gender", "X_intensity", "Y_intensity", "submitted_plate", "submitted_well", "missing_rate", "heterozygosity", "heterozygosity_pc_corrected", "heterozygosity_missing_outlier", "PSCA", "in_kinship", "excluded_kinship_inference", "excess_relatives", "white_british", "pca_calculation", paste0("PC", seq(1,40)), "phasing_autosome", "phasing_X", "phasing_Y")

# Sample eids; This corresponds to a list of all sample eids (with sex information), retrieved from a ".fam" file available from the UKBB portal
sample_eid <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/plink/_001_ukb_cal_chr1_v2.fam", header = F, select = c(1,5), col.names = c("eid", "sex")))


#################################################
### Merge Main Data #############################

df <- cbind(sample_eid, sample_qc)
print(paste0("Start: ", nrow(df), " individuals"))


#################################################
### Sample Filtering ############################

#################################################
### STEP 1: Exclude related samples (pca_calculation = 1)
df <- df[which(df$pca_calculation == 1), ]
print(paste0("STEP 1: Exclude related samples: ", nrow(df), " individuals"))


#################################################
### STEP 2: Exclude non-white, non-British ancestry samples (white_british = 1)
df <- df[which(df$white_british == 1), ]  
print(paste0("STEP 3: Exclude non-white, non-British ancestry samples: ", nrow(df), " individuals"))


#################################################
### STEP 3: Exclude redacted samples & sample for which CNVs were not called (eid = 6026310)
df <- df[-which(df$eid < 0),]  
df <- df[-which(df$eid == 6026310), ]
print(paste0("STEP 4: Exclude redacted samples: ", nrow(df), " individuals"))


#################################################
### STEP 4: Exclude retracted samples; This file was available from the UKBB portal 
retracted <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/org/w16389_20210809.csv", header = F, col.names = "eid"))
df <- df[!df$eid %in% retracted$eid, ]
print(paste0("STEP 5: Exclude retracted samples: ", nrow(df), " individuals"))


################################################
### STEP 5: Exclude genotype plate outliers (samples genotyped on plates with a mean CNV count per sample > 100); Sample eids are contained in "plate_outliers.txt"
plate_outlier <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/raw/plate_outliers.txt", header = F, select = c(1), col.names = "eid"))
df <- df[!df$eid %in% plate_outlier$eid, ]
print(paste0("STEP 6: Exclude genotype plate outliers: ", nrow(df), " individuals"))


#################################################
### STEP 6: Exclude extreme CNV profiles (samples with > 200 CNVs or 1 CNV > 10Mb); "ukb_cnv.gz" is the final PennCNV output, with all CNV calls in a linear format

# CNVs > 10Mb
cnv <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/cnv_calls/PennCNV/final/ukb_cnv.gz", header = F, select = c(3,5), col.names = c("Length", "File")))
cnv$Length <- as.numeric(gsub("length=|,", "", cnv$Length))
length_outlier <- cnv[which(cnv$Length > 10000000), "File"]
print(paste0("There are ", length(unique(length_outlier)), " samples with at least 1 CNV > 10Mb"))

# > 200 CNVs
cnv_sqc <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/cnv_calls/PennCNV/final/ukb_cnv_sqc.gz", header = T, select = c(1,10)))
num_outlier <- cnv_sqc[which(cnv_sqc$NumCNV > 200), "File"]
print(paste0("There are ", length(unique(num_outlier)), " samples with > 200 CNVs"))

# Overlap of both criteria
cnv_outlier <- as.numeric(gsub("_.*", "", union(length_outlier, num_outlier)))
print(paste0("There are ", length(unique(cnv_outlier)), " CNV outlier samples"))

# Exclude individuals
df <- df[!df$eid %in% cnv_outlier, ]
print(paste0("STEP 7: Exclude extreme CNV profiles (>200 CNVs or or 1 CNV >10Mb): ", nrow(df), " individuals"))


#################################################
### STEP 7: Exclude samples with non-matching submitted vs. inferred sex
df <- df[which(df$submitted_gender == df$inferred_gender), ]
print(paste0("STEP 8: Exclude sex mismatches: ", nrow(df), " individuals"))


#################################################
### Save Data ###################################

# All
fwrite(df[, "eid", drop = F], "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/01_samples/data/final/samples_white_british_All.txt", col.names = F, row.names = F, quote = F, sep = "\t")

# Males only
fwrite(df[which(df$sex == 1), "eid", drop = F], "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/01_samples/data/final/samples_white_british_M.txt", col.names = F, row.names = F, quote = F, sep = "\t")

# Females only
fwrite(df[which(df$sex == 2), "eid", drop = F], "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/01_samples/data/final/samples_white_british_F.txt", col.names = F, row.names = F, quote = F, sep = "\t")
