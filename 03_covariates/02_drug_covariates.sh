#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=drug_covariates   	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          					# Define number of nodes required
#SBATCH --time=0-00:10:00   				# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  				# Define the partition on which the job shall run. May be omitted
#SBATCH --mem=15GB           				# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/03_covariates/data/log/02_drug_covariates-%j.out


#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/03_covariates/script/02_drug_covariates.R

echo "Job done"
