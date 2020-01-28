#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --job-name='expected_het'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/expected_het.stdout
#SBATCH -p koeniglab

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
# first standardize all file formats on the command line
# cut PARENTS.txt -f1-2,5-6 > parents.calls
fileNames <- c("parents.calls","F50.calls","F27.calls","F18.calls","F28.calls","F58.calls")
outNames <- c("parents.expected_het","F50.expected_het","F27.expected_het","F18.expected_het","F28.expected_het","F58.expected_het")
# ("F50.expected_het"...)

for (i in 1:length(fileNames)) {
  sample <- read_delim(fileNames[i], delim="\t", col_names = FALSE,, trim_ws = TRUE)
  sample$countsum <- (sample$X3 + sample$X4)
  sample$heterozygosity <- 2*(sample$X3/sample$countsum)*(sample$X4/sample$countsum)
  write.table(sample, file=OutNames[i], quote=F ,sep="\t",row.names=F,col.names=F)
}
