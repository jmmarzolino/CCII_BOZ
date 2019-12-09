#!/bin/bash -l

#SBATCH --ntasks=5
#SBATCH --mem-per-cpu=200G
#SBATCH --time=6-00:00:00
#SBATCH --job-name="gvcf and call"
#SBATCH --output=/rhome/jmarz001/bigdata/CCXXI_POOL/scripts/cmd_003_gvcf.stdout
#SBATCH -p highmem
#SBATCH --array=1-8

# load modules
module load java gatk/4.1.1.0 bcftools/1.8
# load reference genome and the text file with
INDEX=/rhome/jmarz001/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta
CHR=$PROJECT_DIR/args/gatk_chr
INTERVAL=$(head -n $SLURM_ARRAY_TASK_ID $CHR | tail -n 1)
# set directories
PROJECT_DIR=/rhome/jmarz001/bigdata/CCXXI_POOL
SEQS=`$PROJECT_DIR/args/bam_files | cut -f1)`
BAMS=${PROJECT_DIR}/data/bams
cd ${BAMS}
SNPS=${PROJECT_DIR}/data/calls
FILT=${SNPS}/filter

for file in $SEQS
do
gatk HaplotypeCaller -R $INDEX -L ${INTERVAL} -I $BAM/${file} -O $SNPS/${file}_${INTERVAL}.gvcf -ERC GVCF
done

# combine the per-generation+per-chromosome gvcf files into a CCXXI-per-chromosome database for conversion into vcf files
# GenomicsDBImport takes in one or more single-sample GVCFs and imports data over at least one genomics interval, and outputs a directory containing a GenomicsDB datastore with combined multi-sample data.
# GenotypeGVCFs can then read from the created GenomicsDB directly and output the final multi-sample VCF.
cd $SNPS
for file in $SEQS
do
gatk GenomicsDBImport \
  -V ${file}_${INTERVAL}.gvcf \
  -V ${file}_${INTERVAL}.gvcf \
  --genomicsdb-workspace-path $SNPS/${INTERVAL} \
  --intervals ${INTERVAL}
done


GATK=/opt/linux/centos/7.x/x86_64/pkgs/gatk/4.0.12.0/build/libs/gatk-package-4.0.12.0-20-gf9a2e5c-SNAPSHOT-local.jar
for file in $SEQS
do
java -jar $GATK GenotypeGVCFs \
   --max-alternate-alleles 2 \
   --new-qual \
   -L $INTERVAL \
   -R $INDEX \
   -V gendb://$SNPS/${INTERVAL} \
   --output $SNPS/${INTERVAL}_raw_calls.vcf
done

# get the stats you need!
bcftools stats $SNPS/${INTERVAL}_raw_calls.vcf > $FILT/${INTERVAL}.stats






HAPLOTYPE CALLING (ONE FOR EVERY CHROMOSOME)
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --array=1-10%4
#SBATCH --output=gvcf1_5_%a.out
#SBATCH --error=gvcf1_5_%a.err
#SBATCH --job-name="gvcf1_5"
#SBATCH -p koeniglab
cd /rhome/keelyb/bigdata/CCparents/bams/vcfs
sample=`head -n ${SLURM_ARRAY_TASK_ID} bams5 | tail -n1`
module load gatk
module load htslib/1.9
gatk HaplotypeCaller -R /rhome/keelyb/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta \
-I /rhome/keelyb/bigdata/CCparents/bams/rmdup/${sample}.rmdup.bam --emit-ref-confidence GVCF -L chr1H -O ${sample}.chr1H.raw.g.vcf
MAKING A VARIANT DATABASE
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --array=2-7
#SBATCH --job-name="gdbimport"
#SBATCH -p koeniglab
module load gatk
cd /rhome/keelyb/bigdata/CCparents/bams/vcfs/
gatk --java-options -Xmx12g GenomicsDBImport \
-V 45.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 47.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 48.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 50.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 52.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 56.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 59.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 61.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 64.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 69.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 70.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 198.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 199.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 200.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 201.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 202.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 203.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 204.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 205.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 206.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 207.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 208.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 210.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 211.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 212.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 214.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 215.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 216.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 217.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
-V 218.chr${SLURM_ARRAY_TASK_ID}H.raw.g.vcf \
--genomicsdb-workspace-path CCV_chr${SLURM_ARRAY_TASK_ID}H_database -L chr${SLURM_ARRAY_TASK_ID}H
JOINT VARIANT CALLING
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --array=1-7
#SBATCH --output=callsnpsCCII_%a.out
#SBATCH --job-name="callsnps"
#SBATCH -p koeniglab
module load gatk
cd /rhome/keelyb/bigdata/CCparents/bams/vcfs/
gatk GenotypeGVCFs -R /rhome/keelyb/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta \
-V  gendb://CCII_chr${SLURM_ARRAY_TASK_ID}H_database -O CCII_chr${SLURM_ARRAY_TASK_ID}H.vcf
COMBINE VCFS FROM ALL CHROMOSOMES
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --output=combineCCII.out
#SBATCH --job-name="combine"
#SBATCH -p short
module load picard
picard GatherVcfs \
I=CCII_chr1H.vcf \
I=CCII_chr2H.vcf \
I=CCII_chr3H.vcf \
I=CCII_chr4H.vcf \
I=CCII_chr5H.vcf \
I=CCII_chr6H.vcf \
I=CCII_chr7H.vcf \
O=CCII.vcf
