#!/bin/bash -l

#SBATCH --ntasks=5
#SBATCH --mem=200G
#SBATCH --time=6-00:00:00
#SBATCH --job-name="vcf"
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cmd_004_gvcf.stdout
#SBATCH --error=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cmd_004_gvcf.err
#SBATCH -p highmem
#SBATCH --array=1-8

# load modules
module load gatk/4.1.1.0

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
SEQS=$PROJECT_DIR/args/dup_bams
BAMS=${PROJECT_DIR}/data/bams
SNPS=${PROJECT_DIR}/data/calls

# load reference genome and files
INDEX=/rhome/jmarz001/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta
CHR=$PROJECT_DIR/args/gatk_chr
INTERVAL=$(head -n $SLURM_ARRAY_TASK_ID $CHR | tail -n 1)

for file in $SEQS
do
gatk HaplotypeCaller -R $INDEX -I ${BAMS}/${file}.picard_rmdup.bam -L ${INTERVAL} -O $SNPS/${file}_${INTERVAL}.raw.gvcf -ERC GVCF
done
