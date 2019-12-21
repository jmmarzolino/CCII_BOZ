#!/bin/bash -l

#SBATCH --nodes=15
#SBATCH --mem=200G
#SBATCH --job-name="vcf"
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cmd_004_gvcf.stdout
#SBATCH -p highmem
#SBATCH --array=1-8

# load modules
module load samtools/1.9 gatk/4.1.1.0

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
SEQS=$PROJECT_DIR/args/dup_bams
BAMS=${PROJECT_DIR}/data/bams
SNPS=${PROJECT_DIR}/data/calls

# load reference genome and files
INDEX=/rhome/jmarz001/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta
CHR=$PROJECT_DIR/args/gatk_chr
INTERVAL=$(head -n $SLURM_ARRAY_TASK_ID $CHR | tail -n 1)

for file in $(cat < $SEQS)
do
samtools index -c ${BAMS}/${file}.picard_rmdup.bam

gatk HaplotypeCaller -R $INDEX -I ${BAMS}/${file}.picard_rmdup.bam -L ${INTERVAL} -O $SNPS/${file}_${INTERVAL}.raw.gvcf -ERC GVCF
done
