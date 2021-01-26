#!/bin/bash -l
#SBATCH --mem=5G
#SBATCH --time=02:00:00
#SBATCH --output=combine_windows
#SBATCH --job-name="combine_windows"
#SBATCH --partition=short
#

cd ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/POPGENWIN/CALCED_100KB
cat Px* | sort -k1,1 -k2,2n > P.txt
cat F18x* | sort -k1,1 -k2,2n > F18.txt
cat F28x* | sort -k1,1 -k2,2n > F28.txt
cat F58x* | sort -k1,1 -k2,2n > F58.txt

cd ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/POPGENWIN/CALCED_500KB
cat Px* | sort -k1,1 -k2,2n > P.txt
cat F18x* | sort -k1,1 -k2,2n > F18.txt
cat F28x* | sort -k1,1 -k2,2n > F28.txt
cat F58x* | sort -k1,1 -k2,2n > F58.txt

cd ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/POPGENWIN/CALCED_1MB
cat Px* | sort -k1,1 -k2,2n > P.txt
cat F18x* | sort -k1,1 -k2,2n > F18.txt
cat F28x* | sort -k1,1 -k2,2n > F28.txt
cat F58x* | sort -k1,1 -k2,2n > F58.txt

cd ~/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/POPGENWIN/CALCED_5MB
cat Px* | sort -k1,1 -k2,2n > P.txt
cat F18x* | sort -k1,1 -k2,2n > F18.txt
cat F28x* | sort -k1,1 -k2,2n > F28.txt
cat F58x* | sort -k1,1 -k2,2n > F58.txt
