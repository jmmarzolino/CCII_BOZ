#!/usr/bin/env Rscript

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/data/bams/call_variants/")
library(readr)

#F50 <- read_table2("F50.calls", col_names = FALSE)
#F27 <- read_table2("F27.calls", col_names = FALSE)
#F18 <- read_table2("F18.calls", col_names = FALSE)
#F28 <- read_table2("F28.calls", col_names = FALSE)
#F58 <- read_table2("F58.calls", col_names = FALSE)

fileNames <- c("F50.calls","F27.calls","F18.calls","F28.calls","F58.calls")
fileNamesOut <- sub(".calls", "_het.calls" , fileNames)
# fileNamesOut
# ("F50_het.calls","F27_het.calls","F18_het.calls","F28_het.calls","F58_het.calls")

for (i in fileNames) {
  sample <- read_table2(i, col_names = FALSE)
  sample$heterozygosity <- 2* (sample$X3/ (sample$X3+ sample$X4)) *  (sample$X4/ (sample$X3+ sample$X4))
  sample$p <- sample$X3 / (sample$X3+sample$X4)
  sample$q <- 1 - sample$p
}
  write.table(sample, file=fileNamesOut[i], quote=F ,sep="\t",row.names=F,col.names=F)
