#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --job-name='count_filter'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/count_filter.stdout
#SBATCH -p koeniglab
#SBATCH --time=5-00:00:00

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
progeny_counts <- read_delim("progeny_counts","\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)

# remove lines with too high count numbers (sum of row > 1000)
# ie. don't copy lines that are out of bounds
# column numbers are 7 through 16

#find the sum of the row, make it the 17th column
# if the value in the 17th column is less than 1000, copy the row (columns 1 to 16) to the new data.frame
progeny_counts[,17] <- rowSums(progeny_counts[,7:16], na.rm=TRUE)
filtered_counts <- progeny_counts[which(progeny_counts[,17] < 1000),1:16]

write.table(filtered_counts, file="filtered_counts", quote=F ,sep="\t",row.names=F,col.names=F)
