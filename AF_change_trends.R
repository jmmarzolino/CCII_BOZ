#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --time=3-00:00:00
#SBATCH --job-name='fishers'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/AFchange.stdout
#SBATCH -p batch

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
all_freqs_sub <- read_delim("all_freqs_sub",
    "\t", escape_double = FALSE, col_names = FALSE,
    trim_ws = TRUE)

# first, try to burn out ggplot by plotting the dotplots for every genome position


# subset the data to include only the position & allele frequencies from th top one percent of p-values


# get allele frequencies from the subset



library(readr)
all_freqs_sub <- read_delim("all_freqs_sub", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

all_freqs_sub$X9 <- all_freqs_sub$X4

colnames(all_freqs_sub) <- c("chr","position","parents","F18","F27", "F28", "F50", "F58", "F18")

parent_freq
F18_freq
F28_freq
F58_freq

F27_freq
F50_freq

F18inc <- 0
F18dec <- 0
# increase or decrease from parental frequency?
for (i in 1:1000) {
  parent_freq <- all_freqs_sub[i,3]
  F18_freq <- all_freqs_sub[i,4]
  if (F18_freq > parent_freq) {
    F18inc=F18inc+1
  } else {
    F18dec=F18dec+1
  }
}
print("number of F18 sites that increased")
print(F18inc)
print("number of F18 sites that decreased")
print(F18dec)
# is mean a better idea?



if ( test_expression1) {
statement1
} else if ( test_expression2) {
statement2
} else if ( test_expression3) {
statement3
} else {
statement4
}
