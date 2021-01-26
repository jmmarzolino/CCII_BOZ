#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --job-name='cumulative'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cum_pos_freqs.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
all_freqs <- read_delim("all_freqs",
    "\t", escape_double = FALSE, col_names = FALSE,
    trim_ws = TRUE)

for (row in 1:length(all_freqs$X1)) {
  chromosomes <- c("chr1H","chr2H","chr3H","chr4H","chr5H","chr6H","chr7H")
  chr_val <- gsub("chr(\\w+)H", "\\1", all_freqs[row,1])
  for (i in chromosomes[1:7]){
    len1 <- sum(all_freqs[,1]==chromosomes[1])
    len2 <- sum(all_freqs[,1]==chromosomes[2]) + len1
    len3 <- sum(all_freqs[,1]==chromosomes[3]) + len2
    len4 <- sum(all_freqs[,1]==chromosomes[4]) + len3
    len5 <- sum(all_freqs[,1]==chromosomes[5]) + len4
    len6 <- sum(all_freqs[,1]==chromosomes[6]) + len5
    len7 <- sum(all_freqs[,1]==chromosomes[7]) + len6

    if (chr_val=="1"){
      all_freqs[row,4] <- all_freqs[row,9]
    }
    if (chr_val=="2"){
      all_freqs[row,4] <- all_freqs[row,9] + len1
    }
    if (chr_val=="3"){
      all_freqs[row,4] <- all_freqs[row,9] + len2
    }
    if (chr_val=="4"){
      all_freqs[row,4] <- all_freqs[row,9] + len3
    }
    if (chr_val=="5"){
      all_freqs[row,4] <- all_freqs[row,9] + len4
    }
    if (chr_val=="6"){
      all_freqs[row,4] <- all_freqs[row,9] + len5
    }
    if (chr_val=="7"){
      all_freqs[row,4] <- all_freqs[row,9] + len6
    }}}

write.table(all_freqs, file="freqs_and_pos", quote=F ,sep="\t",row.names=F,col.names=F)

