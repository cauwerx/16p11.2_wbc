#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=16p11.2_WBC   			# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          					# Define number of nodes required
#SBATCH --time=0-00:10:00   				# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  				# Define the partition on which the job shall run. May be omitted
#SBATCH --mem=2GB           				# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/05_associations/data/log/02_drug_individual_association-%j.out


#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/05_associations/script/02_drug_individual_association.R

echo "Job done"
