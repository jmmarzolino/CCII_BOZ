#!/bin/bash -l

#SBATCH --ntasks=8
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --job-name='grep'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/grep.out
#SBATCH -p koeniglab

zcat /rhome/jmarz001/shared/SEQ_RUNS/10_8_2018/FASTQ/250_S229_L003_R1_001.fastq.gz | wc -l
echo "250 L003 raw"
zcat /rhome/jmarz001/bigdata/CCII_BOZ/data/trim_fqs/250_L003_1.fq.gz | wc -l
echo "250 L003 trim"

zcat /rhome/jmarz001/shared/SEQ_RUNS/10_8_2018/FASTQ/255_S230_L003_R1_001.fastq.gz | wc -l
echo "255 L003 raw"
zcat /rhome/jmarz001/bigdata/CCII_BOZ/data/trim_fqs/255_L003_1.fq.gz | wc -l
echo "255 L003 trim"
