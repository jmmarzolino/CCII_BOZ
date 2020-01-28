#!/bin/bash -l
#SBATCH --mem=5G
#SBATCH --time=02:00:00
#SBATCH --array=1-87%87
#SBATCH --output=popwins_100kb
#SBATCH --job-name="popwins_100kb"
#SBATCH --partition=short
#
module load python
module load samtools
#
cd ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/POPGENWIN/SUBBED_500KB
mkdir ../CALCED_500KB

thisfile=`ls | head -n ${SLURM_ARRAY_TASK_ID} | tail -n1`
for i in `cat $thisfile`
do
python ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/SCRIPTS/popgenwin.py <(~/bigdata/LOCAL_SOFTWARE/htslib/htslib/tabix ../../CALL_VARIANTS/FULL_AFS.txt.gz $i | cut -f1,2,3,4,5,6) $i >> ../CALCED_500KB/P${thisfile}
python ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/SCRIPTS/popgenwin.py <(~/bigdata/LOCAL_SOFTWARE/htslib/htslib/tabix ../../CALL_VARIANTS/FULL_AFS.txt.gz $i | cut -f1,2,3,4,7,8) $i >> ../CALCED_500KB/F18${thisfile}
python ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/SCRIPTS/popgenwin.py <(~/bigdata/LOCAL_SOFTWARE/htslib/htslib/tabix ../../CALL_VARIANTS/FULL_AFS.txt.gz $i | cut -f1,2,3,4,9,10) $i >> ../CALCED_500KB/F28${thisfile}
python ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/SCRIPTS/popgenwin.py <(~/bigdata/LOCAL_SOFTWARE/htslib/htslib/tabix ../../CALL_VARIANTS/FULL_AFS.txt.gz $i | cut -f1,2,3,4,11,12) $i >> ../CALCED_500KB/F58${thisfile}
done
