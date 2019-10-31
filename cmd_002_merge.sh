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

cd ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/BAM
outfile=`head -n ${SLURM_ARRAY_TASK_ID} ../../INPUT/mergelist.txt | cut -d" " -f1  | tail -n1`
infile1=`head -n ${SLURM_ARRAY_TASK_ID} ../../INPUT/mergelist.txt | cut -d" " -f2  | tail -n1`
infile2=`head -n ${SLURM_ARRAY_TASK_ID} ../../INPUT/mergelist.txt | cut -d" " -f3  | tail -n1`
#infile3=`head -n ${SLURM_ARRAY_TASK_ID} ../../INPUT/mergelist.txt | cut -d" " -f4  | tail -n1`
#infile4=`head -n ${SLURM_ARRAY_TASK_ID} ../../INPUT/mergelist.txt | cut -d" " -f5  | tail -n1`

samtools merge -@10 $outfile $infile1 $infile2 $infile3 $infile4
samtools index -c $outfile
