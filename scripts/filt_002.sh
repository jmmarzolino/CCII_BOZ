#!/bin/bash -l

#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=20G
#SBATCH --time=80:00:00
#SBATCH --job-name='filt snps'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/filt_001.stdout
#SBATCH -p koeniglab
#SBATCH --array=1-4

## PROJECT
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
# load modules
module load bcftools/1.8

# set directories
SNPS=$PROJECT_DIR/data/calls
FILT=$SNPS/filter
mkdir $FILT
# define file list
SEQS=$PROJECT_DIR/args/vcf_files
cd $SNPS ; ls *.vcf > $SEQS
# Define files to run over
FILE=$(head -n $SLURM_ARRAY_TASK_ID $SEQS | tail -n 1)
# ${sample_name}_rawsnps.vcf -> ${sample_name}
sample_name=$(basename "$FILE" | cut -d_ -f2-3)

# bcftools filter [options] <file>
# INCLUDE: sites missing less than 50%, depth higher than 0 and less than 4031,
# N_ALT=1 makes only options ref/single alt allele (biallelic site instead of 4 alt alleles or something)
bcftools view -i 'F_MISSING<0.5 && DP>0 & DP<4031 && QUAL>30 && N_ALT=1' $SNPS/$FILE -o $FILT/${sample_name}_filt001.vcf
# && COUNT(GT="het")<124 && COUNT(GT="alt")>0 && COUNT(GT="ref")>0

