#!/bin/bash -l

#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=5-00:00:00
#SBATCH --job-name="Read Test"
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/BamRead.stdout
#SBATCH -p koeniglab

module load samtools
cd /rhome/jmarz001/bigdata/CCII_BOZ/data/bams

echo "Is read from 250 L003 in merged bam?"
samtools view 250_merged.bam | grep "A00351:60:HCCNTDSXX:3:1523:23818:2785"
echo "Is read from 250 L004 in merged bam?"
samtools view 250_merged.bam | grep "A00351:60:HCCNTDSXX:4:2278:18566:4679"


echo "Is read from 255 L003 in merged bam?"
samtools view 255_merged.bam | grep "A00351:60:HCCNTDSXX:3:1214:23439:28651"
echo "Is read from 255 L004 in merged bam?"
samtools view 255_merged.bam | grep "A00351:60:HCCNTDSXX:4:2472:2492:35556"

cd /rhome/jmarz001/bigdata/CCII_BOZ/data/trim_fqs
echo "Are reads from 250 L003 bam in 250 L003 FQ?"
zcat 250_L003_1.fq.gz | grep "A00351:60:HCCNTDSXX:3:1523:23818:2785"
zcat 250_L003_2.fq.gz | grep "A00351:60:HCCNTDSXX:3:1523:23818:2785"

echo "Are reads from 250 L004 bam in 250 L004 FQ?"
zcat 250_L004_1.fq.gz | grep "A00351:60:HCCNTDSXX:4:2278:18566:4679"
zcat 250_L004_2.fq.gz | grep "A00351:60:HCCNTDSXX:4:2278:18566:4679"

echo "Are reads from 255 L003 bam in 255 L003 FQ 1 or 2?"
zcat 255_L003_1.fq.gz | grep "A00351:60:HCCNTDSXX:3:1214:23439:28651"
zcat 255_L003_2.fq.gz | grep "A00351:60:HCCNTDSXX:3:1214:23439:28651"

echo "Are reads from 255 L004 bam in 255 L004 FQ 1 or 2?"
zcat 255_L004_1.fq.gz | grep "A00351:60:HCCNTDSXX:4:2472:2492:35556"
zcat 255_L004_2.fq.gz | grep "A00351:60:HCCNTDSXX:4:2472:2492:35556"
