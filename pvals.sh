#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --job-name="cut pvals"
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/pvals.stdout
#SBATCH -p batch

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
cd ${PROJECT_DIR}/results

cut -f1-2,17 Genome_Counts_Fishers_values* > pvals
