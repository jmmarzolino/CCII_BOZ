#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --time=02:00:00
#SBATCH --job-name='MAF'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/MAF_plots.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(ggplot2)
options(stringsAsFactors = F)

MAF <- read_delim("MAF_all","\t", col_names = FALSE, trim_ws = TRUE)

MAF_narm <- MAF[which(!is.na(MAF)),]

ggplot(MAF_narm, aes(X3)) + geom_histogram()+ theme_minimal()+xlab("parent minor AFS")
ggplot(MAF_narm, aes(X4)) + geom_histogram()+ theme_minimal()+xlab("F18 AFS")
ggplot(MAF_narm, aes(X5)) + geom_histogram()+ theme_minimal()+xlab("F27 AFS")
ggplot(MAF_narm, aes(X6)) + geom_histogram()+ theme_minimal()+xlab("F28 AFS")
ggplot(MAF_narm, aes(X7)) + geom_histogram()+ theme_minimal()+xlab("F50 AFS")
ggplot(MAF_narm, aes(X8)) + geom_histogram()+ theme_minimal()+xlab("F58 AFS")
