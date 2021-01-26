library(readr)
library(ggplot2)
setwd("/rhome/jmarz001/bigdata/CCII_BOZ/results")

## Import CC II BOZ II minor allele frequencies
BOZ_II_MAF <- read_delim("minor_frequencies_BOZ_II", "\t", col_names = T)
## Import CC II Dav + Boz I minor allele frequencies
BOZ_I_MAF <- read_delim("freqs_cum_pos", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)


## Inner join all the allele frequencies by cumulative count column
tmp <- inner_join(BOZ_I_MAF, BOZ_II_MAF, by=CUM)




## Correlations
cor(F27_II$p, F50$p, method="spearman", use = "pairwise.complete.obs")



