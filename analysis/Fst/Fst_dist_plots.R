#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --time=02:00:00
#SBATCH --job-name='FST_dist'
#SBATCH --output=Fst_dist_plots.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(ggplot2)
options(stringsAsFactors = F)

df <- read_delim("Fst","\t", col_names = T, trim_ws = TRUE)
list <- c("F18-to-parent","F27-to-parent","F28-to-parent","F50-to-parent","F58-to-parent")
binlist <- c(5,10,15, 20, 30, 40, 50)

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
