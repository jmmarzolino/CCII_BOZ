#!/bin/bash -l

#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --time=9-00:00:00
#SBATCH --job-name="recall"
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/parental_site_recall.stdout
#SBATCH -p koeniglab
#SBATCH --array=1-6

# load modules
module load bcftools/1.9
module load samtools/1.9
source activate pyenv

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
SEQS=${PROJECT_DIR}/args/sorted_bams.txt
BAMS=${PROJECT_DIR}/data/bams
RESULTS=${PROJECT_DIR}/results

# how the sites were determined in the first place
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ~/bigdata/BARLEY_CCII_PARENTS_WGS/DATA/OUTPUT/FINALCALLS/FILTER_SNPSONLY.vcf > ../../INPUT/RECALL_FILT_NAMES.bed

# now use those sites to extract allele frequencies in the progeny
# refer to .bed file containing sites
SITES_FILE=/rhome/dkoenig/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/INPUT/RECALL_FILT_NAMES.bed

# get filenames from list
FILE=`head -n ${SLURM_ARRAY_TASK_ID} $SEQS | tail -n1`
NAME=`basename $FILE | cut -d_ -f1`

python $PROJECT_DIR/scripts/extractsite_counts.py $BAMS/$FILE $SITES_FILE > $RESULTS/${NAME}.calls

#cd ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/BAM
#~/bigdata/BARLEY_CCII_POOLSEQ_WGS/SCRIPTS/001_ALLELE_COUNT.py ~/bigdata/BARLEY_CCII_PARENTS_WGS/DATA/OUTPUT/FINALCALLS/FILTER_SNPSONLY.vcf > ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/CALL_VARIANTS/PARENTS.txt
