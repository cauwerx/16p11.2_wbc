# 16p11.2 BP4-BP5 CNV association with white blood cell phenotypes
Code repository for *"Possible association of 16p11.2 copy number variation with altered lymphocyte and neutrophil counts"*

**Contact:** Chiara Auwerx (chiara.auwerx -at- unil.ch) or Giuliana Giannuzzi (giuliana.giannuzzi -at- unimi.it).


## Description of content: 

**01_samples:** Contains the scripts to select unrelated, white British individuals used as the base population for this analysis. 


**02_CNV_carriers:** Contains the scripts to detect 16p11.2 BP4-BP5 CNV carriers in the base population and retrieve data on their drug usage. 


**03_covariates:** Contains the scripts to extract covariates, including classical covariates (sex, age, age^2, batch, PC1-40) and drug usage covariates.


**04_phenotypes:** Contains the scripts to extract phenotypic data (monocyte count, lymphocyte count, WBC count, eosinophil count, neutrophil count), inverse normal transform the phenotypes, and correct them for classical covariates. 


**05_associations:** Contains the scripts to perform association analyses between the selected white blood cell-related phenotypes and the 16p11.2 BP4-BP5 CNV status, according to several CNV-effect models and in a sex-agnostic and sex-specific way. Associations are either done correcting solely for classical covariates or correcting for both classical covariates and drug usage, as a potential confounder of the association.  
