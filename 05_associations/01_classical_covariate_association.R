# Associations between 16p11.2 and WBC using classical covariates

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)


#################################################
### STEP 1: Load data ###########################

# Genotypes; File is generated from 01_samples/01_sample_filtering.R
geno <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/data/final/16p11.2_CNV_profile_HC.txt"))

# Phenotypes INT and corrected by sex; Files from 04_phenotypes/01_phenotype_extraction.R 
pheno_all <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/INT_age_age2_batch_PCs/WBC_WB_INT_age_age2_sex_batch_PCs_All.txt"))
pheno_M <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/INT_age_age2_batch_PCs/WBC_WB_INT_age_age2_batch_PCs_M.txt"))
pheno_F <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/final/INT_age_age2_batch_PCs/WBC_WB_INT_age_age2_batch_PCs_F.txt"))


#################################################
### STEP 2: Create sex-specific dataframes ######

# All individuals
df_all <- left_join(geno, pheno_all, by = "eid") 

# Males only
df_M <- left_join(geno[which(geno$sex == "M"), ], pheno_M, by = "eid")

# Female only
df_F <- left_join(geno[which(geno$sex == "F"), ], pheno_F, by = "eid")

# Combine them in a list
data <- list(df_all, df_M, df_F)


#################################################
### STEP 3: Association studies #################

# Create an empty dataframe to store results of the linear regressions
lm_results <- data.frame()

# Counter
i <- 1

# Loop over sexes
for (sex in 1:3) {

	# Define the dataframe to work with
	df <- data[[sex]]

	# Loop over models
	for (m in c("DUP", "DEL", "MIRROR", "U_SHAPE")) {

		# Loop over phenotypes
		for (p in c("monocyte_count", "lymphocyte_count", "WBC_count", "eosinophil_count", "neutrophil_count")) {

			# Fit the linear regression		
			fit <- lm(df[, p] ~ df[, m], na.action = na.exclude)

			# Fill in the result table
			if (sex == 1) {lm_results[i, "SEX"] <- "All"} else if (sex == 2) {lm_results[i, "SEX"] <- "M"} else if (sex == 3) {lm_results[i, "SEX"] <- "F"} 
			lm_results[i, "MODEL"] <- m
			lm_results[i, "PHENO"] <- p
			lm_results[i, "BETA"] <- summary(fit)$coefficients[2,1]
			lm_results[i, "SE"] <- summary(fit)$coefficients[2,2]
			lm_results[i, "t"] <- summary(fit)$coefficients[2,3]
			lm_results[i, "P"] <- summary(fit)$coefficients[2,4]

			# Increase the counter
			i <- i + 1
		}
	}
}

print(paste0("Nominally significant associations: ", nrow(lm_results[which(lm_results$P <= 0.05), ])))
print(paste0("Significant associations (corrected for number of phenotypes): ", nrow(lm_results[which(lm_results$P <= 0.05/10), ])))


#################################################
### STEP 4: Test for sex differences ############
  
# Create an empty dataframe to store results
sex_differences <- data.frame()

# Counter
i <- 1

# Loop over models
for (m in c("DUP", "DEL", "MIRROR", "U_SHAPE")) {

	# Loop over phenotypes
	for (p in c("monocyte_count", "lymphocyte_count", "WBC_count", "eosinophil_count", "neutrophil_count")) {

		# Check if there is a trend towards a significant effect in the sex combined analysis; if so, test for significant difference across sexes
		if (lm_results[which(lm_results$SEX == "All" & lm_results$MODEL == m & lm_results$PHENO == p), "P"] <= 0.5) {
		
			# Calculate t-statistic for difference in association strength between sexes and associated p-value
			t_stat <- (lm_results[which(lm_results$SEX == "M" & lm_results$MODEL == m & lm_results$PHENO == p), "BETA"] - lm_results[which(lm_results$SEX == "F" & lm_results$MODEL == m & lm_results$PHENO == p), "BETA"])/sqrt((lm_results[which(lm_results$SEX == "M" & lm_results$MODEL == m & lm_results$PHENO == p), "SE"])^2 + (lm_results[which(lm_results$SEX == "F" & lm_results$MODEL == m & lm_results$PHENO == p), "SE"])^2)
			p_val <- 2*pnorm(-abs(t_stat), mean = 0, sd = 1)

			# Fill in the result table
			sex_differences[i, "MODEL"] <- m
			sex_differences[i, "PHENO"] <- p
			sex_differences[i, "t"] <- t_stat
			sex_differences[i, "P"] <- p_val

		# If no trend is present, fill the table with NAs
		} else {

			# Fill in the result table
			sex_differences[i, "MODEL"] <- m
			sex_differences[i, "PHENO"] <- p
			sex_differences[i, "t"] <- NA
			sex_differences[i, "P"] <- NA
		}

		# Increase the counter
		i <- i + 1		
	}
}

print(paste0("Nominally significant sex differences in effect: ", nrow(sex_differences[which(sex_differences$P <= 0.05 & !is.na(sex_differences$P)), ])))
mtc <- 0.05/ nrow(na.omit(sex_differences))
print(paste0("Significant associations (corrected for number of phenotypes + models; p < ", round(mtc, 5),"): ", nrow(sex_differences[which(sex_differences$P <= mtc & !is.na(sex_differences$P)), ])))


#################################################
### Save ########################################

# Linear regression results
fwrite(lm_results, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/05_associations/data/final/classical_covariates/lm_results_16p11.2_WBC_classical_covariates.txt", col.names = T, row.names = F, quote = F, sep = "\t")

# Sex differences
fwrite(sex_differences, "/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/05_associations/data/final/classical_covariates/sex_differences_16p11.2_WBC_classical_covariates.txt", col.names = T, row.names = F, quote = F, sep = "\t")
