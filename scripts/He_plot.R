#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --time=02:00:00
#SBATCH --job-name='He_plot'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/results/He_plot.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(ggplot2)


all <- read_delim("all.expected_het", "\t",col_names = FALSE,trim_ws = TRUE)

generation <- c("F1", "F18", "F27", "F28","F50", "F58")
binlist <- c(5,10,15, 20, 30, 40, 50)

j=1
xlab <- paste(generation[j], "Expected Heterozygosity", sep=" ")
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_Het_expect",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(all, aes(X5)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2750000)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=2
xlab <- paste(generation[j], "Expected Heterozygosity", sep=" ")
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_Het_expect",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(all, aes(X6)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2750000)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=3
xlab <- paste(generation[j], "Expected Heterozygosity", sep=" ")
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_Het_expect",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(all, aes(X7)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2750000)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=4
xlab <- paste(generation[j], "Expected Heterozygosity", sep=" ")
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_Het_expect",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(all, aes(X8)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2750000)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=5
xlab <- paste(generation[j], "Expected Heterozygosity", sep=" ")
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_Het_expect",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(all, aes(X9)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2750000)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=6
xlab <- paste(generation[j], "Expected Heterozygosity", sep=" ")
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_Het_expect",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(all, aes(X10)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2750000)
  ggsave(outname, width = 10, height = 8, units = "in")
}

