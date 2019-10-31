#!/bin/bash -l
#SBATCH --mem-per-cpu=1G
#SBATCH --time=10:00:00
#SBATCH --array=1-2%2
#SBATCH --output=RECALL.stdout
#SBATCH --job-name="RECALL"
#SBATCH --partition=koeniglab
module load bcftools
module load samtools
#
cd ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/BAM
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ~/bigdata/BARLEY_CCII_PARENTS_WGS/DATA/OUTPUT/FINALCALLS/FILTER_SNPSONLY.vcf > ../../INPUT/RECALL_FILT_NAMES.bed
#
thisfileA=`head -n ${SLURM_ARRAY_TASK_ID} ../../INPUT/BAMFILES2MAP.txt | cut -f1 | tail -n1`
thisfileB=`echo $thisfileA | tr '/' '\n' | tail -n1`
#samtools index -c $thisfileA
#mkdir ../CALL_VARIANTS/
python ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/SCRIPTS/extractsite_counts.py $thisfileA ../../INPUT/RECALL_FILT_NAMES.bed > ../CALL_VARIANTS/${thisfileB}.calls



cd ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/BAM
~/bigdata/BARLEY_CCII_POOLSEQ_WGS/SCRIPTS/001_ALLELE_COUNT.py ~/bigdata/BARLEY_CCII_PARENTS_WGS/DATA/OUTPUT/FINALCALLS/FILTER_SNPSONLY.vcf > ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/CALL_VARIANTS/PARENTS.txt
