#!/usr/bin/env Rscript

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
genome_counts <- read_delim("progeny_counts",
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

for (row in 1:nrow(merge)) {
  # take allele counts from each line in turn, convert into a matrix with generations as rows
  site_counts <- matrix(as.numeric(merge[row,3:12]), nrow=5,byrow = TRUE)
  r=12
  for (i in 1:ncol(parent_comp)) {
    # make significance comparisons between the parental counts and each generation in turn
    p_value <- fisher.test(site_counts[i,1:2],alternative = "two.sided")[1] # conduct fishers on whichever row (generation) and get ref and allele counts (1:2)
    merge[row,r] <- c(p_value) #assign the fishers p-value to be added to the end of the current row in order (r)
    r=r+1 
  }
  for (i in pairwise_comp){

  }
  for (i in all6pops){

  }
  for (i in allprogeny){

  }
}

write.table(genome_counts, file="Genome_Counts_Fishers_values", quote=F ,sep="\t",row.names=F,col.names=F)
