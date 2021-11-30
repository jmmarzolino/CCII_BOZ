#!/bin/bash -l

#SBATCH --ntasks=2
#SBATCH --mem=30G
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/BOZ_II/cmd_003_dups1-2.stdout
#SBATCH --job-name="dups"
#SBATCH --time=3-00:00:00
#SBATCH -p koeniglab
#SBATCH --array=1-2

module load samtools/1.9
INDEX=/rhome/jmarz001/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
BAMS=${PROJECT_DIR}/data/bams
SEQS=${PROJECT_DIR}/args/coverage_bams.txt

# get filenames from list
FILE=$(head -n $SLURM_ARRAY_TASK_ID $SEQS | tail -n 1 | cut -f1)
sample_name=$(head -n ${SLURM_ARRAY_TASK_ID} $SEQS | tail -n 1 | cut -f2)

samtools view -b -T $INDEX $FILE | samtools sort -@ 10 > $BAMS/${sample_name}_sorted.bam
samtools index -c $BAMS/${sample_name}_sorted.bam

mkdir -pv /scratch/jmarz/
java -Xmx6g -jar /rhome/jmarz001/picard/build/libs/picard-2.21.3-SNAPSHOT-all.jar MarkDuplicates I=$BAMS/${sample_name}_sorted.bam O=$BAMS/${sample_name}.picard_rmdup.bam M=$BAMS/${sample_name}.metrics.txt REMOVE_DUPLICATES=true TMP_DIR=/scratch/jmarz/
rmdir -r /scratch/jmarz/
