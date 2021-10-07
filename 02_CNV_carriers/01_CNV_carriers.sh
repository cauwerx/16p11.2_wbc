#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=cnv_carriers   	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-00:10:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which the job shall run. May be omitted
#SBATCH --mem=1GB           		# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/data/log/01_CNV_carriers-%j.out


#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/COLLABORATION/Giuliana/16p11.2_wbc/02_CNV_carriers/script/01_CNV_carriers.R

echo "Job done"
