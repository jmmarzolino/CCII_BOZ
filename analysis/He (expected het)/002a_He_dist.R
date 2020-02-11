#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --time=02:00:00
#SBATCH --job-name='He_plot'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/He_plot.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(ggplot2)
# load in the He and weighted_He files
all_He <- read_delim("all_He", "\t", col_names=T,trim_ws = TRUE)
all_weighted_He <- read_delim("all_weighted_He", "\t", col_names=T,trim_ws = TRUE)
gens <- c("F0","F18","F27","F28","F50","F58")

# generate distribution plots for
for (m in gens) {
  OutName = paste(paste(m,"He",sep="_"),"dist.pdf",sep="_")
  OutName2 = paste(paste(m,"weighted_He",sep="_"),"dist.pdf",sep="_")

  xlab <- paste(m, "Expected Heterozygosity", sep=" ")
  xlab2 <- paste(m, "Weighted Expected Heterozygosity", sep=" ")

  #plotting He
  ggplot(all_He, aes(get(m))) + geom_histogram(bins=15)+ theme_minimal()+xlab(xlab)+ylim(0,2750000)+xlim(0,0.5)
  ggsave(OutName, width = 10, height = 8, units = "in")

  #plotting weighted He
  ggplot(all_weighted_He, aes(get(m))) + geom_histogram(bins=15)+ theme_minimal()+xlab(xlab2)+xlim(0,0.5)
  ggsave(OutName2, width = 10, height = 8, units = "in")
}
