#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --output=cmd_006_combine.stdout
#SBATCH --job-name="combine"
#SBATCH -p short

module load picard/2.18.3
cd $SNPS

picard GatherVcfs \
I=chr1H_raw_calls.vcf \
I=chr2H_raw_calls.vcf \
I=chr3H_raw_calls.vcf \
I=chr4H_raw_calls.vcf \
I=chr5H_raw_calls.vcf \
I=chr6H_raw_calls.vcf \
I=chr7H_raw_calls.vcf \
I=chrUn_raw_calls.vcf \
O=CCII_DavBoz.vcf
