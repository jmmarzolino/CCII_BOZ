#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --time=02:00:00
#SBATCH --job-name='FST_dist'
#SBATCH --output=/bigdata/koeniglab/jmarz001/CCII_BOZ/scripts/Fst_dist_plots.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(ggplot2)
options(stringsAsFactors = F)

df <- read_delim("Fst_cumpos","\t", col_names = T, trim_ws = TRUE)
df <- df[,3:7]
list <- colnames(df)
binlist <- c(10,15,20,25)

j=1
xlab <- list(paste(list[j], "FST", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(list[j], "FST",sep="_"),m,sep="_"),"bins.png",sep="")
  col <- list[j]
  ggplot(df, aes(col)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=2
xlab <- list(paste(list[j], "FST", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(list[j], "FST",sep="_"),m,sep="_"),"bins.png",sep="")
  col <- list[j]
  ggplot(df, aes(col)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=3
xlab <- list(paste(list[j], "FST", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(list[j], "FST",sep="_"),m,sep="_"),"bins.png",sep="")
  col <- list[j]
  ggplot(df, aes(col)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=4
xlab <- list(paste(list[j], "FST", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(list[j], "FST",sep="_"),m,sep="_"),"bins.png",sep="")
  col <- list[j]
  ggplot(df, aes(col)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=5
xlab <- list(paste(list[j], "FST", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(list[j], "FST",sep="_"),m,sep="_"),"bins.png",sep="")
  col <- list[j]
  ggplot(df, aes(col)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)
  ggsave(outname, width = 10, height = 8, units = "in")
}

