#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --time=02:00:00
#SBATCH --job-name='het_plot'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/het_plot.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(ggplot2)

F18_het <- read_delim("F18_het.calls","\t", col_names = FALSE, trim_ws = TRUE)
F28_het <- read_delim("F28_het.calls","\t", col_names = FALSE, trim_ws = TRUE)
F58_het <- read_delim("F58_het.calls","\t", col_names = FALSE, trim_ws = TRUE)
F27_het <- read_delim("F27_het.calls","\t", col_names = FALSE, trim_ws = TRUE)
F50_het <- read_delim("F50_het.calls","\t", col_names = FALSE, trim_ws = TRUE)
Fparent_het <- read_delim("parent_het.calls","\t", col_names = FALSE, trim_ws = TRUE)

F18_narm <- F18_het[which(!is.na(F18_het$X7)),]
F27_narm <- F27_het[which(!is.na(F27_het$X7)),]
F28_narm <- F28_het[which(!is.na(F28_het$X7)),]
F50_narm <- F50_het[which(!is.na(F50_het$X7)),]
F58_narm <- F58_het[which(!is.na(F58_het$X7)),]
parent_narm <- Fparent_het[which(!is.na(Fparent_het$X9)),]

ggplot(F18_narm, aes(X7)) + geom_histogram()+ theme_minimal()+xlab("F18 Expected Heterozygosity")
ggplot(F27_narm, aes(X7)) + geom_histogram()+ theme_minimal()+xlab("F27 Expected Heterozygosity")
ggplot(F28_narm, aes(X7)) + geom_histogram()+ theme_minimal()+xlab("F28 Expected Heterozygosity")
ggplot(F50_narm, aes(X7)) + geom_histogram()+ theme_minimal()+xlab("F50 Expected Heterozygosity")
ggplot(F58_narm, aes(X7)) + geom_histogram()+ theme_minimal()+xlab("F58 Expected Heterozygosity")
ggplot(parent_narm, aes(X9)) + geom_histogram()+ theme_minimal()+xlab("parents Expected Heterozygosity")
