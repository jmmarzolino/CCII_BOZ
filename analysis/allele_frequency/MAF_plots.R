#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --time=02:00:00
#SBATCH --job-name='MAF'
#SBATCH --output=MAF_plots.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(ggplot2)
options(stringsAsFactors = F)

MAF <- read_delim("minor_frequencies","\t", col_names = FALSE, trim_ws = TRUE)

MAF_narm <- MAF[which(!is.na(MAF)),]

ggplot(MAF_narm, aes(X3)) + geom_histogram()+ theme_minimal()+xlab("F1 minor AFS")
ggplot(MAF_narm, aes(X4)) + geom_histogram()+ theme_minimal()+xlab("F18 AFS")
ggplot(MAF_narm, aes(X5)) + geom_histogram()+ theme_minimal()+xlab("F27 AFS")
ggplot(MAF_narm, aes(X6)) + geom_histogram()+ theme_minimal()+xlab("F28 AFS")
ggplot(MAF_narm, aes(X7)) + geom_histogram()+ theme_minimal()+xlab("F50 AFS")
ggplot(MAF_narm, aes(X8)) + geom_histogram()+ theme_minimal()+xlab("F58 AFS")


#########################

generation <- c("F1", "F18", "F27","F28", "F50","F58")
binlist <- c(5,10,15, 20, 30, 40, 50)

j=1
xlab <- list(paste(generation[j], "Allele Frequency Spectrum", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(MAF, aes(X3)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=2
xlab <- list(paste(generation[j], "Allele Frequency Spectrum", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(MAF, aes(X4)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=3
xlab <- list(paste(generation[j], "Allele Frequency Spectrum", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(MAF, aes(X5)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=4
xlab <- list(paste(generation[j], "Allele Frequency Spectrum", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(MAF, aes(X6)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=5
xlab <- list(paste(generation[j], "Allele Frequency Spectrum", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(MAF, aes(X7)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)
  ggsave(outname, width = 10, height = 8, units = "in")
}

j=6
xlab <- list(paste(generation[j], "Allele Frequency Spectrum", sep=" "))
for (m in binlist) {
  outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
  ggplot(MAF, aes(X8)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)
  ggsave(outname, width = 10, height = 8, units = "in")
}
