#!/bin/bash -l

#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=20G
#SBATCH --time=168:00:00
#SBATCH --job-name='call snps'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cmd_002_call.stdout
#SBATCH -p koeniglab
#SBATCH --array=1-4

## PROJECT
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ

# load modules
module load freebayes/1.2.0
# define variable locations
INDEX=/rhome/jmarz001/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta

# set directories
BAMS=$PROJECT_DIR/data/bams
SNPS=$PROJECT_DIR/data/calls
FILT=$SNPS/filter
# create data/file directory structure
mkdir $SNPS
mkdir $FILT

SEQS=$PROJECT_DIR/args/bam_files
cd $BAMS ; ls *.bam > $SEQS
# Define files to run over
FILE=$(head -n $SLURM_ARRAY_TASK_ID $SEQS | tail -n 1)
# 250_L003.bam -> 267_188
sample_name=$(basename "$FILE" | cut -d. -f1)

#freebayes -f [reference] [infiles.bam] > [outfiles.vcf]
freebayes -k -f $INDEX $BAMS/$FILE > $SNPS/${sample_name}_rawsnps.vcf

# get the stats you need!
bcftools stats $SNPS/${sample_name}_rawsnps.vcf > $FILT/${sample_name}.stats
