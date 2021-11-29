#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --job-name='cumulative'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cum_pos.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
pvals <- read_delim("pvals",
    "\t", escape_double = FALSE, col_names = FALSE,
    trim_ws = TRUE)

for (row in 1:length(pvals$X1)) {
  chromosomes <- c("chr1H","chr2H","chr3H","chr4H","chr5H","chr6H","chr7H")
  chr_val <- gsub("chr(\\w+)H", "\\1", pvals[row,1])
  for (i in chromosomes[1:7]){
    len1 <- sum(pvals[,1]==chromosomes[1])
    len2 <- sum(pvals[,1]==chromosomes[2]) + len1
    len3 <- sum(pvals[,1]==chromosomes[3]) + len2
    len4 <- sum(pvals[,1]==chromosomes[4]) + len3
    len5 <- sum(pvals[,1]==chromosomes[5]) + len4
    len6 <- sum(pvals[,1]==chromosomes[6]) + len5
    len7 <- sum(pvals[,1]==chromosomes[7]) + len6

    if (chr_val=="1"){
      pvals[row,4] <- pvals[row,2]
    }
    if (chr_val=="2"){
      pvals[row,4] <- pvals[row,2] + len1
    }
    if (chr_val=="3"){
      pvals[row,4] <- pvals[row,2] + len2
    }
    if (chr_val=="4"){
      pvals[row,4] <- pvals[row,2] + len3
    }
    if (chr_val=="5"){
      pvals[row,4] <- pvals[row,2] + len4
    }
    if (chr_val=="6"){
      pvals[row,4] <- pvals[row,2] + len5
    }
    if (chr_val=="7"){
      pvals[row,4] <- pvals[row,2] + len6
    }}}

write.table(pvals, file="pval_and_pos", quote=F ,sep="\t",row.names=F,col.names=F)
