#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --time=1-00:00:00
#SBATCH --job-name='combine_counts'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/AF000_combine_counts.stdout
#SBATCH -p koeniglab

# set up the environment
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)

# read in the parent call file to set the chr+pos columns
all_counts <- read_delim("F0.calls","\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)

# list every file that needs to be read in
file_list <- c("F18.calls", "F27.calls", "F28.calls", "F50.calls", "F58.calls")
# loop to add each file in turn (by generation NUMERIC ORDER) to the combined file
for (file in file_list){
  Bobby <- read_delim(file,"\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
  all_counts <- cbind(all_counts[1:10000,],Bobby[1:10000,3:4])
}

# add column names for clarity
colnames(all_counts) <- c("CHR","POS","F0 REF", "F0 ALT","F18 REF","F18 ALT","F27 REF","F27 ALT","F28 REF","F28 ALT","F50 REF","F50 ALT","F58 REF","F58 ALT")
# write out the combined count file
write.table(all_counts, file="all_counts", quote=F ,sep="\t",row.names=F,col.names=T)

