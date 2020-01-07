#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --job-name="gvcfUn"
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cmd_004_gvcfUn.stdout
#SBATCH -p highmem
#SBATCH --array=1-5
#SBATCH --time=9-00:00:00

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
INTERVAL=$(tail -n1 $CHR)
FILE=$(head -n $SLURM_ARRAY_TASK_ID $SEQS | tail -n 1)

samtools index -c ${BAMS}/${FILE}.picard_rmdup.bam

gatk HaplotypeCaller -R $INDEX -I ${BAMS}/${FILE}.picard_rmdup.bam -L ${INTERVAL} -O $SNPS/${FILE}_${INTERVAL}.raw.gvcf -ERC GVCF
