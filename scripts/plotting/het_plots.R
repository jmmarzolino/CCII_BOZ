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
# load file
files <- c("F18_het.calls", "F28_het.calls", "F58_het.calls", "F27_het.calls", "F50_het.calls", "parent_het.calls")
names <- c("F18_het", "F28_het", "F58_het", "F27_het", "F50_het", "parent_het")

for (i in 1:6){
  df<- read_delim(files[i],"\t", col_names = FALSE, trim_ws = TRUE)
  df_naremove <- df[which(!is.na(df$X7)),]

  ggplot(df_naremove, aes(X7)) + geom_histogram(binwidth=5)+ theme_minimal()

  # +xlab("genome position")+ylab("p-value")
  # save the plot
  outname <- paste(names[i], ".pdf", sep="")
  ggsave(outname)
}

