# Extract most recent phenotype data, INT, and correct for classical covariates

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)


#################################################
### STEP 1: Load Main Data ######################

# This corresponds to a .csv file with two columns: field ID and short name for the corresponding phenotype
# Phenotypes: monocyte count, lymphocyte count, WBC count, eosinophil count, neutrophil count
pheno_list <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/raw/WBC_traits.csv", header = T))


#################################################
### STEP 2: Order raw phenotype files ###########

# As phenotypic data is spread throughout multiple files downloaded from the UKBB portal, we start by listing all phenotype files available and order them by date
raw_pheno_file <- file.info(list.files("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno", pattern = "^ukb.*csv", full.names = T, recursive = F))
raw_pheno_file <- data.frame(File = rownames(raw_pheno_file), Date = raw_pheno_file[,4])
raw_pheno_file <- raw_pheno_file[order(raw_pheno_file$Date, decreasing = T), ]


#################################################
### STEP 3: Identify most recent location #######
# As mentioned in STEP 2, there are multiple phenotype files, with some phenotypes being in several files
# We use phenotype information originating from the most recent phenotype file

pheno_list$File <- NA
pheno_list$Col <- NA

# Loop over all phenotypes to detect the most recent location
for (p in 1:nrow(pheno_list)) {

	# Define the phenotype
	pheno <- pheno_list[p, "Pheno"] 
	ID <- pheno_list[p, "FieldID"] 
	print(paste0("Extracting most recent file for ", pheno))
	counter <- 1

	# Determine the most recent location --> loop through raw files
	while (is.na(pheno_list[p, "File"])) {
	
		header <- fread(as.character(raw_pheno_file[counter, "File"]), nrow = 0)
		col <- grep(paste0("^", ID, "-"), names(header))
		if (length(col) > 0) {
			pheno_list[p, "File"] <- as.character(raw_pheno_file[counter, "File"])
			print(paste0("Most recent location: ", sub(".*/", "", pheno_list[p, "File"])))
			pheno_list[p, "Col"] <- paste(col, collapse = "_")}
		counter <- counter +1
	}
} 
rm(p, pheno, ID, counter, header, col)
fwrite(pheno_list, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/pheno_list.txt", col.names = T, row.names = F, quote = F, sep = "\t")


#################################################
### STEP 4: Extract phenotypes and average ######

# This file contains all UKBB sample eids, except for the redacted ones
phenotypes <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/general_data/UKBB/eid_no_redacted3.txt.gz", header = T, select = c(2)), drop = F)

# Loop over all phenotypes to extract the data and calculate the average over different measurement instances
for(p in 1:nrow(pheno_list)) {
	
	# Define the phenotype
	pheno <- pheno_list[p, "Pheno"] 
	ID <- pheno_list[p, "FieldID"] 
	file <- pheno_list[p, "File"] 
	col <- pheno_list[p, "Col"] 
	print(paste0("Extracting ", pheno, " from ", sub(".*/", "", file)))

	# Extract columns corresponding to the defined phenotype
	temp_pheno <- as.data.frame(fread(file, header = T, select = c(1, as.numeric(str_split(col, pattern = "_")[[1]]))))

	# Average (if more than one column present)
	if (ncol(temp_pheno) > 2) {
	temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = rowMeans(temp_pheno[, 2:ncol(temp_pheno)], na.rm = T))}

	# Merge
	colnames(temp_pheno) <- c("eid", pheno)
	phenotypes <- full_join(phenotypes, temp_pheno, by = "eid")

}
rm(pheno, ID, file, col, temp_pheno)
print(paste0("Dimensions of phenotype table: ", ncol(phenotypes), " pheno x ", nrow(phenotypes), " eids"))


#################################################
### STEP 5: Select samples ######################
# Retain unrelated, white British samples selected with 01_samples/01_sample_filtering.R

# All individuals
eid_all <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/01_samples/data/final/samples_white_british_All.txt", header = F, col.names = "eid"))
pheno_all <- phenotypes[phenotypes$eid %in% eid_all$eid, ]
print(paste0("Dimensions of phenotype table for all unrelated white British individuals: ", ncol(pheno_all), " pheno x ", nrow(pheno_all), " eids"))
fwrite(pheno_all, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/raw/WBC_WB_raw_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

# Males only
eid_M <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/01_samples/data/final/samples_white_british_M.txt", header = F, col.names = "eid"))
pheno_M <- phenotypes[phenotypes$eid %in% eid_M$eid, ]
print(paste0("Dimensions of phenotype table for unrelated white British males: ", ncol(pheno_M), " pheno x ", nrow(pheno_M), " males"))
fwrite(pheno_M, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/raw/WBC_WB_raw_M.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

# Females only
eid_F <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/01_samples/data/final/samples_white_british_F.txt", header = F, col.names = "eid"))
pheno_F <- phenotypes[phenotypes$eid %in% eid_F$eid, ]
print(paste0("Dimensions of phenotype table for unrelated white British females: ", ncol(pheno_F), " pheno x ", nrow(pheno_F), " females"))
fwrite(pheno_F, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/raw/WBC_WB_raw_F.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


#################################################
### STEP 6: Summary #############################

# All individuals
summary_pheno_all <- data.frame(Pheno = colnames(pheno_all)[2:ncol(pheno_all)],
								Num = colSums(!is.na(pheno_all[, c(2:ncol(pheno_all))])),
								NumNA = colSums(is.na(pheno_all[, c(2:ncol(pheno_all))])),
								Mean = colMeans(pheno_all[, c(2:ncol(pheno_all))], na.rm = T),
								Median = colMedians(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T),
								SD = colSds(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T),
								Min = colMins(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T),
								Max = colMaxs(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T))
fwrite(summary_pheno_all, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/raw/summary_WBC_WB_raw_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

# Males only
summary_pheno_M <- data.frame(Pheno = colnames(pheno_M)[2:ncol(pheno_M)],
							  Num = colSums(!is.na(pheno_M[, c(2:ncol(pheno_M))])),
							  NumNA = colSums(is.na(pheno_M[, c(2:ncol(pheno_M))])),
							  Mean = colMeans(pheno_M[, c(2:ncol(pheno_M))], na.rm = T),
							  Median = colMedians(as.matrix(pheno_M[, c(2:ncol(pheno_M))]), na.rm = T),
							  SD = colSds(as.matrix(pheno_M[, c(2:ncol(pheno_M))]), na.rm = T),
							  Min = colMins(as.matrix(pheno_M[, c(2:ncol(pheno_M))]), na.rm = T),
							  Max = colMaxs(as.matrix(pheno_M[, c(2:ncol(pheno_M))]), na.rm = T))
fwrite(summary_pheno_M, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/raw/summary_WBC_WB_raw_M.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

# Females only
summary_pheno_F <- data.frame(Pheno = colnames(pheno_F)[2:ncol(pheno_F)],
							  Num = colSums(!is.na(pheno_F[, c(2:ncol(pheno_F))])),
							  NumNA = colSums(is.na(pheno_F[, c(2:ncol(pheno_F))])),
							  Mean = colMeans(pheno_F[, c(2:ncol(pheno_F))], na.rm = T),
							  Median = colMedians(as.matrix(pheno_F[, c(2:ncol(pheno_F))]), na.rm = T),
							  SD = colSds(as.matrix(pheno_F[, c(2:ncol(pheno_F))]), na.rm = T),
							  Min = colMins(as.matrix(pheno_F[, c(2:ncol(pheno_F))]), na.rm = T),
							  Max = colMaxs(as.matrix(pheno_F[, c(2:ncol(pheno_F))]), na.rm = T))
fwrite(summary_pheno_F, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/raw/summary_WBC_WB_raw_F.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


#################################################
### STEP 7: Inverse Normal Transformation #######

# All individuals
pheno_int_all <- data.frame(eid = pheno_all$eid)
for (p in 2:ncol(pheno_all)) {
	pheno <- colnames(pheno_all)[p]
	print(paste0("INT - All: ", pheno))
	pheno_int_all[, pheno] <- qnorm((rank(pheno_all[, p], na.last = "keep")-0.5)/sum(!is.na(pheno_all[, p])))
}; rm (p, pheno)
fwrite(pheno_int_all, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/INT/WBC_WB_INT_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

# Males only
pheno_int_M <- data.frame(eid = pheno_M$eid)
for (p in 2:ncol(pheno_M)) {
	pheno <- colnames(pheno_M)[p]
	print(paste0("INT - Males: ", pheno))
	pheno_int_M[, pheno] <- qnorm((rank(pheno_M[, p], na.last = "keep")-0.5)/sum(!is.na(pheno_M[, p])))
}; rm (p, pheno)
fwrite(pheno_int_M, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/INT/WBC_WB_INT_M.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

# Females only
pheno_int_F <- data.frame(eid = pheno_F$eid)
for (p in 2:ncol(pheno_F)) {
	pheno <- colnames(pheno_F)[p]
	print(paste0("INT - Females: ", pheno))
	pheno_int_F[, pheno] <- qnorm((rank(pheno_F[, p], na.last = "keep")-0.5)/sum(!is.na(pheno_F[, p])))
}; rm (p, pheno)
fwrite(pheno_int_F, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/INT/WBC_WB_INT_F.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


#################################################
### STEP 8: Correct for covariates ##############

# Read in "classical" covariates from 03_covariates/01_classical_covariates.R
cov_all <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/03_covariates/data/final/classical_covariates.txt", header = T)) # 335954 eids

# Select subsets (sex is not used when assessing sex-specific traits)
cov_M <- cov_all[cov_all$eid %in% eid_M$eid, -c(4)] 
cov_F <- cov_all[cov_all$eid %in% eid_F$eid, -c(4)]

# Regress out covariates - All individuals
pheno_int_cor_all <- data.frame(eid = pheno_int_all$eid)
for (p in 2:ncol(pheno_int_all)) {
	# Define the phenotype
	pheno <- colnames(pheno_int_all)[p]
	print(paste0("Covariates - All: ", pheno))
	# Create a temporary dataframe
	temp <- right_join(cov_all, pheno_int_all[, c(1,p)], by = "eid")
	colnames(temp)[ncol(temp)] <- "pheno"
	# Regress out covariates
	pheno_int_cor_all[, pheno] <- residuals(lm(pheno ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)
fwrite(pheno_int_cor_all, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/INT_age_age2_batch_PCs/WBC_WB_INT_age_age2_sex_batch_PCs_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

# Regress out covariates - Males only
pheno_int_cor_M <- data.frame(eid = pheno_int_M$eid)
for (p in 2:ncol(pheno_int_M)) {
	# Define the phenotype
	pheno <- colnames(pheno_int_M)[p]
	print(paste0("Covariates - Males: ", pheno))
	# Create a temporary dataframe
	temp <- right_join(cov_M, pheno_int_M[, c(1,p)], by = "eid")
	colnames(temp)[ncol(temp)] <- "pheno"
	# Regress out covariates
	pheno_int_cor_M[, pheno] <- residuals(lm(pheno ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)
fwrite(pheno_int_cor_M, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/INT_age_age2_batch_PCs/WBC_WB_INT_age_age2_batch_PCs_M.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

# Regress out covariates - Females only
pheno_int_cor_F <- data.frame(eid = pheno_int_F$eid)
for (p in 2:ncol(pheno_int_F)) {
	# Define the phenotype
	pheno <- colnames(pheno_int_F)[p]
	print(paste0("Covariates - Females: ", pheno))
	# Create a temporary dataframe
	temp <- right_join(cov_F, pheno_int_F[, c(1,p)], by = "eid")
	colnames(temp)[ncol(temp)] <- "pheno"
	# Regress out covariates
	pheno_int_cor_F[, pheno] <- residuals(lm(pheno ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)
fwrite(pheno_int_cor_F, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/INT_age_age2_batch_PCs/WBC_WB_INT_age_age2_batch_PCs_F.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
