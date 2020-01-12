#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --job-name='fishers'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/fishers_priority.stdout
#SBATCH -p koeniglab
#SBATCH --time=5-00:00:00

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(matrixStats)
progeny_counts <- read_delim("progeny_counts",
    "\t", escape_double = FALSE, col_names = FALSE,
    trim_ws = TRUE)

# remove lines with too high count numbers ((>/< than mean +/- 4SDs)
# first define the upper and lower bounds of the counts for the row
# next, loop through every row, checking for counts beyond the defined bounds
# if ALL of the counts for this row are within bounds (less than upper & more than lower), then copy that line to the next line of the filtered df
# ie. don't copy lines that are out of bounds
# column numbers are 7 through 16
newrow <- 1
storage <- c(0)
progeny_counts_filtered <- c(1:16)
dim(progeny_counts_filtered) <- c(1,16)
progeny_counts_filtered <- data.frame(progeny_counts_filtered)

#if (all(progeny_counts[1,7:ncol(progeny_counts)] < upper & progeny_counts[1,7:ncol(progeny_counts)] > lower)) {progeny_counts_filtered[1,] <- progeny_counts[1,]}
#if (all(progeny_counts[2,7:ncol(progeny_counts)] < upper & progeny_counts[2,7:ncol(progeny_counts)] > lower)) {progeny_counts_filtered[2,] <- progeny_counts[2,]}

for (row in 1:nrow(progeny_counts)){
  storage <- c(as.numeric(progeny_counts[1,7]),as.numeric(progeny_counts[1,8]),as.numeric(progeny_counts[1,9]),as.numeric(progeny_counts[1,10]),as.numeric(progeny_counts[1,11]),as.numeric(progeny_counts[1,12]),as.numeric(progeny_counts[1,13]),as.numeric(progeny_counts[1,14]),as.numeric(progeny_counts[1,15]),as.numeric(progeny_counts[1,16]))
  upper <- mean(storage) + (4*sd(storage))
  lower <- mean(storage) - (4*sd(storage))
  if (all(progeny_counts[row,7:ncol(progeny_counts)] < upper & progeny_counts[row,7:ncol(progeny_counts)] > lower)) {
    progeny_counts_filtered[newrow,] <- progeny_counts[row,]
    newrow <- newrow+1
  }
}
write.table(progeny_counts_filtered, file=progeny_counts_filtered, quote=F ,sep="\t",row.names=F,col.names=F)

# define rows in the matrix to grab for pairwise (and bigger) comparisons
parents <- 1
F18 <- 2
F27 <- 3
F28 <- 4
F50 <- 5
F58 <- 6

#compare everything
all6pops <- c(parents, F18, F27, F28, F50, F58)

# now compute the p-values
for (j in 1:3570){
  print(j)
  end <- j*1000
  start <- end-999
  for (row in start:end) {
    # take allele counts from each line in turn, convert into a matrix with generations as rows
    site_counts <- matrix(as.numeric(progeny_counts[row,5:16]),nrow=6,byrow=T)

    p_value <- fisher.test(site_counts[all6pops,1:2],workspace = 2e8, hybrid=T)[1]
    progeny_counts[row,17] <- p_value
  }
  file_name <- paste("Genome_Counts_Fishers_values_allpops", sep="", j)
  write.table(progeny_counts[start:end,], file=file_name, quote=F ,sep="\t",row.names=F,col.names=F)
}
