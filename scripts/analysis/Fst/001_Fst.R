#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=10:00:00
#SBATCH --job-name='Fst'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/Fst.stdout
#SBATCH -p koeniglab

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
options(stringsAsFactors = F)

df <- read_delim("minor_frequencies","\t", col_names = FALSE, trim_ws = TRUE)

Fst <- df[,1:2]
colnames(Fst) <- c("chr", "pos")

prog <- data.frame(df$X3,df$X4)
colnames(prog) <- c("p_parent", "p_progeny")
prog$p_bar <- (prog$p_parent + prog$p_progeny)/2
prog$Ht <- 2*(prog$p_bar)*(1-prog$p_bar)
prog$Hs <- 2*(prog$p_progeny)*(1-prog$p_progeny)
prog$Fst <- (prog$Ht - prog$Hs)/prog$Ht
Fst$F18 <- abs(prog$Fst)

prog <- data.frame(df$X3,df$X5)
colnames(prog) <- c("p_parent", "p_progeny")
prog$p_bar <- (prog$p_parent + prog$p_progeny)/2
prog$Ht <- 2*(prog$p_bar)*(1-prog$p_bar)
prog$Hs <- 2*(prog$p_progeny)*(1-prog$p_progeny)
prog$Fst <- (prog$Ht - prog$Hs)/prog$Ht
Fst$F27 <- abs(prog$Fst)

prog <- data.frame(df$X3,df$X6)
colnames(prog) <- c("p_parent", "p_progeny")
prog$p_bar <- (prog$p_parent + prog$p_progeny)/2
prog$Ht <- 2*(prog$p_bar)*(1-prog$p_bar)
prog$Hs <- 2*(prog$p_progeny)*(1-prog$p_progeny)
prog$Fst <- (prog$Ht - prog$Hs)/prog$Ht
Fst$F28 <- abs(prog$Fst)

prog <- data.frame(df$X3,df$X7)
colnames(prog) <- c("p_parent", "p_progeny")
prog$p_bar <- (prog$p_parent + prog$p_progeny)/2
prog$Ht <- 2*(prog$p_bar)*(1-prog$p_bar)
prog$Hs <- 2*(prog$p_progeny)*(1-prog$p_progeny)
prog$Fst <- (prog$Ht - prog$Hs)/prog$Ht
Fst$F50 <- abs(prog$Fst)

prog <- data.frame(df$X3,df$X8)
colnames(prog) <- c("p_parent", "p_progeny")
prog$p_bar <- (prog$p_parent + prog$p_progeny)/2
prog$Ht <- 2*(prog$p_bar)*(1-prog$p_bar)
prog$Hs <- 2*(prog$p_progeny)*(1-prog$p_progeny)
prog$Fst <- (prog$Ht - prog$Hs)/prog$Ht
Fst$F58 <- abs(prog$Fst)

write.table(Fst, file="Fst", quote=F ,sep="\t",row.names=F,col.names=T)
