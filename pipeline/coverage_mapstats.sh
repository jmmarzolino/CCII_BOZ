#!/bin/bash -l

#SBATCH --ntasks=8
#SBATCH --mem=20G
#SBATCH --time=168:00:00
#SBATCH --job-name='trim+align'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cmd_001_trim_align.stdout
#SBATCH -p koeniglab
#SBATCH --array=1-5

# load modules
module load samtools/1.9

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
SEQS=${PROJECT_DIR}/args/coverage_bams.txt
BAMS=${PROJECT_DIR}/data/bams

# Define files to run over
FILE=$(head -n $SLURM_ARRAY_TASK_ID $SEQS | tail -n 1 | cut -f1)
generation=$(head -n $SLURM_ARRAY_TASK_ID $SEQS | tail -n 1 | cut -f2)

# Get mapping stats
samtools flagstat $FILE > $BAMS/mappingstats/${generation}_mapstats.txt
