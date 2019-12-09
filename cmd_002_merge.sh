#!/bin/bash -l

#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00:00
#SBATCH --job-name="merge bams"
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cmd_002_merge.stdout
#SBATCH -p batch
#SBATCH --array=1-2

# load modules
module load samtools/1.9
#INDEX=/rhome/jmarz001/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
SEQS=${PROJECT_DIR}/args/merge_bams
BAMS=${PROJECT_DIR}/data/bams
cd $BAMS

# get filenames from list
MERGED_FILE=`head -n ${SLURM_ARRAY_TASK_ID} $SEQS | cut -f3  | tail -n1`
FILE_1=`head -n ${SLURM_ARRAY_TASK_ID} $SEQS | cut -f1  | tail -n1`
FILE_2=`head -n ${SLURM_ARRAY_TASK_ID} $SEQS | cut -f2  | tail -n1`

# merge files and index
samtools merge -@10 $MERGED_FILE $FILE_1 $FILE_2
# would need to parse the merge file name
#samtools view -b -T $INDEX $BAMS/$MERGED_FILE | samtools sort -@ 10 > $BAMS/${MERGED_FILE}_sorted.bam
samtools index -c ${MERGED_FILE}.bam
