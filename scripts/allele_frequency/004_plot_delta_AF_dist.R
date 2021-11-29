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
##############################################################################
#####################Read In Dataset##########################################
##############################################################################
#Load in files with change in allele frequencies calculated
delta_AF_F1toALL <- read_delim("delta_AF_F1toALL","\t", col_names = T, trim_ws = TRUE)
delta_AF_DAVIS <- read_delim("delta_AF_DAVIS","\t", col_names = T, trim_ws = TRUE)
delta_AF_BOZ <- read_delim("delta_AF_BOZ","\t", col_names = T, trim_ws = TRUE)
#bind all the data columns together
df <- cbind.data.frame(delta_AF_F1toALL,delta_AF_DAVIS[,3:5],delta_AF_BOZ[,3:5])
#write_delim(df,"delta_AF_all",delim="\t",col_names=T)
###########################################################################
#Plot actual AF changes as distributions
###########################################################################
df <- read_delim("delta_AF_all","\t",col_names=T,trim_ws=T)
for (x in c(3,8,9,11,12)){
  gen <- colnames(df)[x]
  xlab <- paste("Allele Frequency Change",gen)
  OutName <- paste("delta_AF",gen,sep="_")
  
  g <- ggplot(df, aes(get(gen))) + geom_histogram(bins=40)+theme_minimal() + xlab(xlab)+ylab("Frequency")+ylim(NA,1750000)

  OutName2 <- paste0(OutName, "_distribution.jpeg")
  ggsave(OutName2, g, width=10, height=6, units="in")
}

###########################################################################
#Plot absolute AF changes as distributions
###########################################################################
# edit data frame to be absolute values
#df <- read_delim("delta_AF_all",delim="\t",col_names=T)
df[,3:ncol(df)] <- abs(df[,3:ncol(df)])

for (x in 3:ncol(df)){
  gen <- colnames(df)[x]
  xlab <- paste("Absolute Allele Frequency Change",gen)
  OutName <- paste("delta_AF",gen,sep="_")

  g <- ggplot(df, aes(get(gen))) + geom_histogram(bins=40)+theme_minimal() + xlab(xlab)+ylab("Frequency")

  OutName2 <- paste0(OutName, "_abs_distribution.jpeg")
  ggsave(OutName2, g, width=10, height=6, units="in")
}
