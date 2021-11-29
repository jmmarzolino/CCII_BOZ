#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --job-name='expected_het'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/expected_het.stdout
#SBATCH -p koeniglab
#SBATCH --array=1-6

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# set up environment
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(ggplot2)
options(stringsAsFactors = F)

# list the generations to name out files
gens <- c("F0","F18","F27","F28","F50","F58")
OutName <- paste(gens[task_id],"expected_het",sep="_")
# read in data
df <- read_delim(OutName, delim="\t", col_names = T, trim_ws = T)

# overall contribution of Heterozygosity is weighted by our /certainty/ about genotpye sequencing
binlist <- c(5,10,15, 20, 30, 40, 50)

xlab <- paste(gens[task_id], "Expected Heterozygosity", sep=" ")
for (m in binlist) {
  OutName2 = paste(paste(OutName,m,sep="_"),"bins.pdf",sep="")
  g<-ggplot(df, aes(het)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2750000)
  ggsave(g,OutName2, width = 10, height = 8, units = "in")
}

xlab <- paste(gens[task_id], "Weighted Expected Heterozygosity", sep=" ")
for (m in binlist) {
  OutName2 = paste(paste(paste(OutName,"weighted",sep="_"),m,sep="_"),"bins.pdf",sep="")
  ggplot(df, aes(weighted_het)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2750000)
  ggsave(outname, width = 10, height = 8, units = "in")
}
