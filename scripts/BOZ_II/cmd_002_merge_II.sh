#!/bin/bash -l

#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=50G
#SBATCH --time=1-00:00:00
#SBATCH --job-name="merge bams"
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/BOZ_II/cmd_002_merge.stdout
#SBATCH -p koeniglab
#SBATCH --array=1-2

# load modules
module load samtools/1.9

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
SEQS=${PROJECT_DIR}/args/merge_bams_II
BAMS=${PROJECT_DIR}/data/bams
cd $BAMS

# get filenames from list
MERGED_FILE=`head -n ${SLURM_ARRAY_TASK_ID} $SEQS | cut -f3  | tail -n1`
FILE_1=`head -n ${SLURM_ARRAY_TASK_ID} $SEQS | cut -f1  | tail -n1`
FILE_2=`head -n ${SLURM_ARRAY_TASK_ID} $SEQS | cut -f2  | tail -n1`

# merge files and index
samtools merge -f -@10 $MERGED_FILE $FILE_1 $FILE_2
samtools index -c $MERGED_FILE
