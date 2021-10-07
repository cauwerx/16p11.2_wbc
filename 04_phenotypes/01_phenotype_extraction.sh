#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=pheno_cont   	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-01:00:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which the job shall run. May be omitted
#SBATCH --mem=10GB           		# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/data/log/01_phenotype_extraction-%j.out

#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/04_phenotypes/script/01_phenotype_extraction.R

echo "Job done"
