#!/bin/bash -l

#SBATCH --ntasks=3
#SBATCH --mem-per-cpu=90G
#SBATCH --time=6-00:00:00
#SBATCH --job-name="line count"
#SBATCH --output=lines.stdout
#SBATCH -p intel


module load samtools


samtools index -c F18.picard_rmdup.bam
samtools index -c F28.picard_rmdup.bam
samtools index -c F58.picard_rmdup.bam

samtools view F18.picard_rmdup.bam | wc -l
samtools view F27.picard_rmdup.bam | wc -l
samtools view F28.picard_rmdup.bam | wc -l
samtools view F50.picard_rmdup.bam | wc -l
samtools view F58.picard_rmdup.bam | wc -l
