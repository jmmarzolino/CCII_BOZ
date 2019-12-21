#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --job-name='fishers'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/fishers.stdout
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

#comparisons with parental allele freqs
parent_F18 <- c(parents, F18)
parent_F27 <- c(parents, F27)
parents_F28 <- c(parents, F28)
parents_F50 <- c(parents, F50)
parents_F58 <- c(parents, F58)

# informed pairs
late <- c(F50, F58)
early <- c(F27,F28)
DB_transition <- c(F18, F27)
Davis_earlylate <- c(F18,F58)
Boz_earlylate <- c(F27, F50)

# between sequential generations
F18_F28 <- c(F18, F28)
F28_F58 <- c(F28, F58)
F27_F50 <- c(F27, F50)

parent_comp <- data.frame(parent_F18,parent_F27,parents_F28,parents_F50, parents_F58)
pairwise_comp <- data.frame(late, early, DB_transition, Davis_earlylate, Boz_earlylate, F18_F28, F28_F58, F27_F50)
#compare everything
all6pops <- c(parents, F18, F27, F28, F50, F58)
allprogeny <- c(F18, F27, F28, F50, F58)

for (j in 1:3570){
  print(j)
  end <- j*1000
  start <- end-999
  for (row in start:end) {
    # take allele counts from each line in turn, convert into a matrix with generations as rows
    site_counts <- matrix(as.numeric(progeny_counts[row,5:16]),nrow=6,byrow=T)
    r=12 #number of filled columns in loaded table
    for (i in 2:ncol(parent_comp)) {
      parents <- site_counts[1,1:2]
      progeny <- site_counts[i,1:2]
      comparison <- matrix(c(parents,progeny), nrow=2, byrow=T)
      # make significance comparisons between the parental counts and each generation in turn
      p_value <- fisher.test(comparison,alternative = "two.sided",workspace = 2e8)[1] # conduct fishers on whichever row (generation) and get ref and allele counts (1:2)
      progeny_counts[row,(r+1)] <- p_value #assign the fishers p-value to be added to the end of the current row in order (r)
      r=r+1
    }
    for (i in 1:ncol(pairwise_comp)){
      a <- site_counts[pairwise_comp[1,i],1:2]
      b <- site_counts[pairwise_comp[2,i],1:2]
      comparison <- matrix(c(a,b), nrow=2, byrow=T)
      p_value <- fisher.test(comparison,alternative = "two.sided",workspace = 2e8)[1]
      progeny_counts[row,r] <- p_value
      r=r+1
    }
    p_value <- fisher.test(site_counts[allprogeny,1:2],workspace = 2e8,, hybrid=T)[1]
    progeny_counts[row,r] <- p_value
    r=r+1

    p_value <- fisher.test(site_counts[all6pops,1:2],workspace = 2e8, hybrid=T)[1]
    progeny_counts[row,r] <- p_value
  }
  file_name <- paste("Genome_Counts_Fishers_values", sep="", j)
  write.table(progeny_counts[start:end,], file=file_name, quote=F ,sep="\t",row.names=F,col.names=F)
}
