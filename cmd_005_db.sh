#!/bin/bash -l

#SBATCH --ntasks=5
#SBATCH --mem=200G
#SBATCH --job-name="db"
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cmd_005_db.stdout
#SBATCH -p highmem
#SBATCH --array=1-8

# load modules
module load gatk/4.1.1.0 bcftools/1.8

# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCII_BOZ
SEQS=$PROJECT_DIR/args/dup_bams
SNPS=${PROJECT_DIR}/data/calls
FILT=${SNPS}/filter

# load reference genome and files
INDEX=/rhome/jmarz001/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta
CHR=$PROJECT_DIR/args/gatk_chr
INTERVAL=$(head -n $SLURM_ARRAY_TASK_ID $CHR | tail -n 1)

# combine the per-generation+per-chromosome gvcf files into a per-chromosome database for conversion into vcf files
# GenomicsDBImport takes in one or more single-sample GVCFs and imports data over at least one genomics interval, and outputs a directory containing a GenomicsDB datastore with combined multi-sample data.
# GenotypeGVCFs can then read from the created GenomicsDB directly and output the final multi-sample VCF.

gatk --java-options -Xmx200g GenomicsDBImport \
  -V $SNPS/F18_${INTERVAL}.raw.gvcf \
  -V $SNPS/F27_${INTERVAL}.raw.gvcf \
  -V $SNPS/F28_${INTERVAL}.raw.gvcf \
  -V $SNPS/F50_${INTERVAL}.raw.gvcf \
  -V $SNPS/F58_${INTERVAL}.raw.gvcf \
  --genomicsdb-workspace-path $SNPS/${INTERVAL} \
  --intervals ${INTERVAL}

gatk --java-options -Xmx200g GenotypeGVCFs \
   --max-alternate-alleles 2 \
   --new-qual \
   -L $INTERVAL \
   -R $INDEX \
   -V gendb://$SNPS/${INTERVAL} \
   -O $SNPS/${INTERVAL}_raw_calls.vcf

# get the stats you need!
bcftools stats $SNPS/${INTERVAL}_raw_calls.vcf > $FILT/${INTERVAL}.stats
