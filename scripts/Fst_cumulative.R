#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name='cumulative'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/Fst_cumulative.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(pacman)
p_load(ggplot2, dplyr, tidyr, data.table)
options(stringsAsFactors = F)

df<-fread("Fst")
OutName <- "Fst"

# parse locus
OutName1<-paste(OutName, "cumpos", sep="_")
names(df)[1]<-"CHR"
names(df)[2]<-"POS"

# format for plotting
df$BP<-as.numeric(df$POS)
result <- df %>%
  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(df, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)

write.table(result, file=OutName1, quote=F ,sep="\t",row.names=F,col.names=T)

