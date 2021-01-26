#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --job-name='count_filter'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/AF001_count_filter.stdout
#SBATCH -p koeniglab
#SBATCH --time=50:00:00

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
all_counts <- read_delim("all_counts","\t", escape_double = FALSE, col_names = T,trim_ws = TRUE)

# remove lines with too high count numbers (sum of row > 1000)
# ie. don't copy lines that are out of bounds
# count column numbers are 5 through 14

#find the sum of the row, make it the 15th column
# if the value in the 15th column is less than 1000, copy the row (columns 1 to 16) to the new data.frame
all_counts[,15] <- rowSums(all_counts[,5:14], na.rm=TRUE)
filtered_counts <- all_counts[which(all_counts[,15] < 1000),1:14]

write.table(filtered_counts, file="filtered_counts", quote=F ,sep="\t",row.names=F,col.names=T)
