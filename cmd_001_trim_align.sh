#!/bin/bash -l

#SBATCH --ntasks=8
#SBATCH --mem=20G
#SBATCH --time=168:00:00s
#SBATCH --job-name='trim+align'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cmd_001_trim_align.stdout
#SBATCH -p koeniglab
#SBATCH --array=1-4%4

# general pattern for using trimmomatic on the cluster
###java -jar <path to trimmomatic jar> SE/PE [-threads <threads>] [-phred33/64] [-trimlog <log file>] <input1> <input2> <paired output 1> <unpaired output 1> <paired output 2> <unpaired output 2> <step 1>...

# load modules
module load trimmomatic/0.36
module load minimap2 samtools

# define variable locations
TRIMMOMATIC=/opt/linux/centos/7.x/x86_64/pkgs/trimmomatic/0.36/trimmomatic.jar
ADAPTERDIR=/rhome/cfisc004/software/Trimmomatic-0.36/adapters
INDEX=/rhome/dkoenig/shared/GENOMES/NEW_BARLEY/GENOME_SPLIT/barley_split_reference.fa

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
SEQS=${PROJECT_DIR}/args/raw_fqs
TRIMMED=${PROJECT_DIR}/data/trim_fq
BAMS=${PROJECT_DIR}/data/bams

# create data/file directory structure
mkdir ${PROJECT_DIR}/data
mkdir $TRIMMED
mkdir $BAMS

# get filenames from list
FILE_1=$(head -n $SLURM_ARRAY_TASK_ID $SEQS | tail -n 1 | cut -f1)
FILE_2=$(head -n $SLURM_ARRAY_TASK_ID $SEQS | tail -n 1 | cut -f2)
sample_name=$(head -n ${SLURM_ARRAY_TASK_ID} $FILE_LIST | tail -n1 | cut -f3 | cut -d_ -f1)
run_name=$(head -n ${SLURM_ARRAY_TASK_ID} $FILE_LIST | tail -n1 | cut -f3 | cut -d_ -f3)
# rhome/jmarz001/shared/SEQ_RUNS/10_8_2018/FASTQ/250_S229_L003_R1_001.fastq       rhome/jmarz001/shared/SEQ_RUNS/10_8_2018/FASTQ/250_S229_L003_R2_001.fastq       250_S229_L003_001       ACAGTG  Novaseq

# Quality/Adapter trimming
java -jar $TRIMMOMATIC PE -threads 10 \
$FILE_1 $FILE_2 \
$TRIMMED/${sample_name}_${run_name}_1.fq.gz $TRIMMED/${sample_name}_${run_name}_1.un.fq.gz \
$TRIMMED/${sample_name}_${run_name}_1.fq.gz $TRIMMED/${sample_name}_${run_name}_1.un.fq.gz \
ILLUMINACLIP:"$ADAPTERDIR"/PE_all.fa:2:30:10 \
SLIDINGWINDOW:4:20 MINLEN:75




minimap2 -t 10 -ax sr ~/shared/GENOMES/NEW_BARLEY/2019_Release_Morex_vers2/Barley.mmi -R "@RG\tID:bar_samp_${SLURM_ARRAY_TASK_ID}\tPL:illumina\tSM:${sample_name}\tLB:${sample_name}" ${wkdir}/OUTPUT/TRIM_FASTQS/${sample_name}${run_name}_1.fq.gz ${wkdir}/OUTPUT/TRIM_FASTQS/${sample_name}${run_name}_2.fq.gz > ${wkdir}/OUTPUT/BAM/${sample_name}${run_name}.sam

samtools view -b -T ~/shared/GENOMES/NEW_BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta ${wkdir}/OUTPUT/BAM/${sample_name}${run_name}.sam | samtools sort -@ 20 > ${wkdir}/OUTPUT/BAM/${sample_name}${run_name}.bam
samtools index -c ${wkdir}/OUTPUT/BAM/${sample_name}${run_name}.bam


# mapping stats
mkdir $RES/mappingstats/
samtools flagstat $RES/"$SAMPLE".sam > $RES/mappingstats/"$SAMPLE"_mapstats.txt

# sam to sorted bam and index long bams with csi file
samtools view -bS $RES/"$SAMPLE".sam | samtools sort -o $RES/"$SAMPLE".bam
rm *.sam
samtools index -c $RES/"$SAMPLE".bam

# extract unmapped reads
mkdir $RES/unmapped
samtools view -f4 -b $RES/"$SAMPLE".bam > $RES/"$SAMPLE".unmapped.bam
# export unmapped reads from original reads
bedtools bamtofastq -i $RES/"$SAMPLE".unmapped.bam -fq $RES/"$SAMPLE".unmapped.fq
mv *unmapped* unmapped/
