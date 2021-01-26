#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=02:00:00
#SBATCH --job-name='delta_AF'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/plot_delta_AF_dist.stdout
#SBATCH -p koeniglab

#Set up environment
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(pacman)
p_load(ggplot2, dplyr, tidyr, data.table)
options(stringsAsFactors = F)


df <- read_delim("delta_AF_all","\t",col_names=T,trim_ws=T)
for (x in c(3,8,9,11,12)){
  gen <- colnames(df)[x]
  xlab <- paste("Allele Frequency Change",gen)
  OutName <- paste("delta_AF",gen,sep="_")
  
  g <- ggplot(df, aes(get(gen))) + geom_histogram(bins=40)+theme_minimal() + xlab(xlab)+ylab("Frequency")+ylim(NA,1750000)

  OutName2 <- paste0(OutName, "_distribution.jpeg")
  ggsave(OutName2, g, width=10, height=6, units="in")
}
