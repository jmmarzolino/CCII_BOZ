#!/bin/bash -l

#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/collect_pvals.stdout
#SBATCH --job-name="collect"
#SBATCH --time=3-00:00:00
#SBATCH -p koeniglab

cd /bigdata/koeniglab/jmarz001/CCII_BOZ/results
LENGTH=$(ls -l filtered_Fishers_pvalues_allpops* | wc -l)

cat test{1..$LENGTH} >> pvals
