#!/bin/bash -l

#SBATCH --ntasks=8
#SBATCH --mem=20G
#SBATCH --time=74:00:00
#SBATCH --job-name='trim+align'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cmd_001_trim_align_II.stdout
#SBATCH -p koeniglab
#SBATCH --array=1-4

# general pattern for using trimmomatic
###java -jar <path to trimmomatic jar> SE/PE [-threads <threads>] [-phred33/64] [-trimlog <log file>] <input1> <input2> <paired output 1> <unpaired output 1> <paired output 2> <unpaired output 2> <step 1>...

# load modules
module load trimmomatic/0.36
module load minimap2/2.17 samtools/1.9

# define variable locations
TRIMMOMATIC=/opt/linux/centos/7.x/x86_64/pkgs/trimmomatic/0.36/trimmomatic.jar
ADAPTERDIR=/rhome/cfisc004/software/Trimmomatic-0.36/adapters
INDEX=/rhome/jmarz001/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
RAW_SEQ_DIR=/rhome/jmarz001/shared/SEQ_RUNS/1_5_2021/FASTQ/JMDK01
#ls /rhome/dkoenig/shared/SEQ_RUNS/1_5_2021/FASTQ/JMDK01 > ${PROJECT_DIR}/args/raw_fqs_II #then edited into two columns
SEQS=${PROJECT_DIR}/args/raw_fqs_II
TRIMMED=${PROJECT_DIR}/data/trim_fqs
BAMS=${PROJECT_DIR}/data/bams

# get filenames from list
FILE_1=$(head -n $SLURM_ARRAY_TASK_ID $SEQS | tail -n 1 | cut -f1)
FILE_2=$(head -n $SLURM_ARRAY_TASK_ID $SEQS | tail -n 1 | cut -f2)
sample_name=$(basename $FILE_1 | cut -d_ -f1)

# Quality/Adapter trimming
java -jar $TRIMMOMATIC PE -threads 10 \
$RAW_SEQ_DIR/$FILE_1 $RAW_SEQ_DIR/$FILE_2 \
$TRIMMED/${sample_name}_1.fastq.gz $TRIMMED/${sample_name}_1.un.fastq.gz \
$TRIMMED/${sample_name}_2.fastq.gz $TRIMMED/${sample_name}_2.un.fastq.gz \
ILLUMINACLIP:"$ADAPTERDIR"/PE_all.fa:2:30:10 \
SLIDINGWINDOW:4:20 MINLEN:75

# Map to reference
minimap2 -t 10 -ax sr /rhome/jmarz001/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley.mmi -R "@RG\tID:${sample_name}_${SLURM_ARRAY_TASK_ID}\tPL:illumina\tSM:${sample_name}\tLB:${sample_name}" $TRIMMED/${sample_name}_1.fastq.gz $TRIMMED/${sample_name}_2.fastq.gz > $BAMS/${sample_name}.sam

# Get mapping stats
# this didn't work and I'm not yet sure why
#mkdir $BAMS/mappingstats/
#samtools flagstat $BAMS/${sample_name}.sam > $BAMS/mappingstats/${sample_name}_mapstats.txt

# Convert sam to sorted bam and index bams with csi
samtools view -b -T $INDEX $BAMS/${sample_name}.sam | samtools sort -@ 20 > $BAMS/${sample_name}.bam
samtools index -c $BAMS/${sample_name}.bam
