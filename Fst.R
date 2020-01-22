#!/usr/bin/env Rscript

#SBATCH --ntasks=4
#SBATCH --mem=60G
#SBATCH --time=10-00:00:00
#SBATCH --job-name='Fst'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/Fst.stdout
#SBATCH -p batch

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
options(stringsAsFactors = F)

df <- read_delim("minor_frequencies","\t", col_names = FALSE, trim_ws = TRUE)
outdf[,1:2] <- df[,1:2]
p_parent <- df[,3]

for (i in 4:8){
  p_progeny <- df[,i]
  p_bar <- ((p_parent + p_progeny)/2)
  p_var <- (((p_progeny - p_bar) * (p_progeny - p_bar))/2)
  Fst <- ((p_var)/(p_bar))
  outdf[row,(i-1)] <- Fst
  }

colnames(outdf) <- c("chr", "pos", "F18-to-parent", "F27-to-parent", "F28-to-parent", "F50-to-parent", "F58-to-parent")
write.table(outdf, file="Fst", quote=F ,sep="\t",row.names=F,col.names=T)
