#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --job-name='fishers'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/fishers_priority.stdout
#SBATCH -p koeniglab

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
progeny_counts <- read_delim("progeny_counts",
    "\t", escape_double = FALSE, col_names = FALSE,
    trim_ws = TRUE)

# define rows in the matrix to grab for pairwise (and bigger) comparisons
parents <- 1
F18 <- 2
F27 <- 3
F28 <- 4
F50 <- 5
F58 <- 6

#compare everything
all6pops <- c(parents, F18, F27, F28, F50, F58)

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
