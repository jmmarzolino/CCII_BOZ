#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH --job-name='WindowsAndCount_functions'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/windowed_counts.stdout
#SBATCH -p short

# set up environment
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
options(stringsAsFactors = F)

# to make windows file make a genome file in the format "chr# \t ChrLength" on the command line
#ARGS=/rhome/jmarz001/bigdata/CCII_BOZ/args
#INDEX=/rhome/jmarz001/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta.fai
#awk -v OFS='\t' {'print $1,$2'} $INDEX > $ARGS/barley_bedgenome

# read in allele count file
progeny_counts <- read_delim("filtered_counts","\t", col_names = FALSE, trim_ws = TRUE)




windows <- function()
