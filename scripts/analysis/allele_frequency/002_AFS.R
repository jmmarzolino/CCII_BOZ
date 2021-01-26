#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --time=02:00:00
#SBATCH --job-name='AFS'
#SBATCH --output=/bigdata/koeniglab/jmarz001/CCII_BOZ/results/AFS.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(ggplot2)
options(stringsAsFactors = F)

df <- read_delim("minor_frequencies","\t", col_names = T, trim_ws = TRUE)

generation <- c("F0", "F18", "F27","F28", "F50","F58")
binlist <- c(5,10,15, 20, 30, 40, 50)

# save all the figures to corresponding folder
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results/figures")
j=1
xlab <- list(paste(generation[j], "Allele Frequency Spectrum", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(df, aes(`F0 MAF`)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2500000)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=2
xlab <- list(paste(generation[j], "Allele Frequency Spectrum", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(df, aes(`F18 MAF`)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2500000)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=3
xlab <- list(paste(generation[j], "Allele Frequency Spectrum", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(df, aes(`F27 MAF`)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2500000)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=4
xlab <- list(paste(generation[j], "Allele Frequency Spectrum", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(df, aes(`F28 MAF`)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2500000)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=5
xlab <- list(paste(generation[j], "Allele Frequency Spectrum", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(df, aes(`F50 MAF`)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2500000)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=6
xlab <- list(paste(generation[j], "Allele Frequency Spectrum", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(df, aes(`F58 MAF`)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2500000)
  ggsave(outname, width = 10, height = 8, units = "in")
}
