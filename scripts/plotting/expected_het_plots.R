#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --time=02:00:00
#SBATCH --job-name='het_plot'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/expected_het_plots.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(ggplot2)

# combine the expected_het files into one df which can be plotted by column & various bin no.s
# copy the chr, pos, and He from parent file

fileNames <- c("parents.expected_het","F50.expected_het","F27.expected_het","F18.expected_het","F28.expected_het","F58.expected_het")
all.expected_het <- read_delim(fileNames[1],"\t", col_names = FALSE, trim_ws = TRUE)

# start reading in files starting AFTER parent file [1] which serves as file base
# copy in the progeny file's fifth column (expected het) to a new column
for (i in 2:length(fileNames)) {
  sample <- read_delim(fileNames[i],"\t", col_names = FALSE, trim_ws = TRUE)
  all.expected_het[,ncol(all.expected_het)+1] <- sample[,5]
}
write.table(all.expected_het, file="all.expected_het", quote=F ,sep="\t",row.names=F,col.names=F)

generation <- c("parents", "F50", "F27", "F18", "F28", "F58")
binlist <- c(15, 20, 30, 40, 50)

for (j in 5:ncol(all.expected_het)) {
  xlab <- paste(generation[j], "Expected Heterozygosity", sep=" ")
  for (m in binlist) {
    outname = paste(paste(paste(generation[j], "_AFS",sep=""),m,sep="_"),"bins.png",sep="")
    # generation = F50 + _AFS = F50_AFS
    # + _ + binlist[#]==m  >> F50_AFS_15
    # + "bins.png" >> F50_AFS_15bins.png
    ggplot(all.expected_het, aes(all.expected_het[,j])) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)
    ggsave(outname, width = 10, height = 8, units = "in")
  }
}

