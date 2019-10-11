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
DATA_DIR=${PROJECT_DIR}/data
TRIMMED=${PROJECT_DIR}/data/trim_fq
BAMS=${PROJECT_DIR}/data/bams

# create data/file directory structure
mkdir $DATA_DIR
mkdir $TRIMMED
mkdir $BAMS

# get filenames from list
FILE=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1)
# get basename of file, stripping at "."
#267_S234_L003_R1_001.fastq ==> 267_S234_L003
NAME=$(basename "$FILE" | cut -d. -f1 | cut -d_ -f1-3)
# == 267_L003/4
SHORT=$(basename "$NAME" | cut -d_ -f1,3)

# Quality/Adapter trimming
java -jar $TRIMMOMATIC PE -threads 10 \
$WORKINGDIR/"$NAME"_R1_001.fastq $WORKINGDIR/"$NAME"_R2_001.fastq \
$RESULTSDIR/"$SHORT"_1_trimmed_paired.fq $RESULTSDIR/"$SHORT"_1_unpaired.fq \
$RESULTSDIR/"$SHORT"_2_trimmed_paired.fq $RESULTSDIR/"$SHORT"_2_unpaired.fq \
ILLUMINACLIP:"$ADAPTERDIR"/PE_all.fa:2:30:10 \
SLIDINGWINDOW:4:20 MINLEN:75

${thisfile_1} ${thisfile_2}
${wkdir}/OUTPUT/TRIM_FASTQS/${sample_name}${run_name}_1.fq.gz ${wkdir}/OUTPUT/TRIM_FASTQS/${sample_name}${run_name}_1_un.fq.gz
${wkdir}/OUTPUT/TRIM_FASTQS/${sample_name}${run_name}_2.fq.gz ${wkdir}/OUTPUT/TRIM_FASTQS/${sample_name}${run_name}_2_un.fq.gz


# select each set of reads from the file in turn
thisfile_1=`head -n ${SLURM_ARRAY_TASK_ID} $FILE_LIST | tail -n1 | cut -f1`
thisfile_2=`head -n ${SLURM_ARRAY_TASK_ID} $FILE_LIST | tail -n1 | cut -f2`
sample_name=`head -n ${SLURM_ARRAY_TASK_ID} $FILE_LIST | tail -n1 | cut -f3`
run_name=`head -n ${SLURM_ARRAY_TASK_ID} $FILE_LIST | tail -n1 | cut -f4`


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
