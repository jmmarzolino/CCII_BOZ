#!/bin/bash -l

#SBATCH --ntasks=16
#SBATCH --mem=20G
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cmd_003_dups.stdout
#SBATCH --job-name="dups"
#SBATCH --time=9-00:00:00
#SBATCH -p koeniglab
#SBATCH --array=1-8

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
BAMS=${PROJECT_DIR}/data/bams
SEQS=${PROJECT_DIR}/args/coverage_bams.txt

# get filenames from list
FILE=$(head -n $SLURM_ARRAY_TASK_ID $SEQS | tail -n 1 | cut -f1)
sample_name=$(head -n ${SLURM_ARRAY_TASK_ID} $SEQS | tail -n 1 | cut -f2)

mkdir -pv /scratch/jmarz/picard/
java -Xmx16g -jar /rhome/jmarz001/picard/build/libs/picard-2.21.3-SNAPSHOT-all.jar MarkDuplicates I=${FILE} O=$BAMS/${sample_name}.picard_rmdup.bam M=$BAMS/${sample_name}.metrics.txt REMOVE_DUPLICATES=true TMP_DIR=/scratch/jmarz/picard/
rmdir -r /scratch/jmarz/picard/
