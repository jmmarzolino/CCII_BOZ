#!/bin/bash -l

#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=20G
#SBATCH --time=2:00:00
#SBATCH --job-name='call snps'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/a001_make_windows.stdout
#SBATCH -p short

module load bedtools/2.28.0

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
ARGS=$PROJECT_DIR/args

INDEX=/rhome/jmarz001/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta.fai

# to make windows file, bedtools requires a genome file
# genome file must be in the format "chr# \t LengthOfChrom"
# there can be ONLY two columns, so I will make this file based on the indexed ref genome
awk -v OFS='\t' {'print $1,$2'} $INDEX > $ARGS/barley_bedgenome

bedtools makewindows -g $ARGS/barley_bedgenome -w 100000 > $ARGS/genome.100kbwindows.bed
bedtools makewindows -g $ARGS/barley_bedgenome -w 1000000 > $ARGS/genome.1Mwindows.bed
