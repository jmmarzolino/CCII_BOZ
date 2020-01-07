#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --time=7-00:00:00
#SBATCH --job-name='subset and graph'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/top_1percent.stdout
#SBATCH -p batch

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(ggplot2)

# load table
pval_and_pos <- read_delim("pval_and_pos", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# set the cutoff for p-values at 1%
cutoff <- quantile(pval_and_pos$X3,0.01)
# subset and save the data that is at or below the cutoff
top_1per <- pval_and_pos[which(pval_and_pos$X3 <= cutoff),]
## might need to re-assign cummulative chr order, I'll see after graphing if it works

#vconvert p-vales into positive values and make 1's into 0's with a -log transform!
top_1per$X5 <- -log10(top_1per$X3)
# save the subsetted table by writing out
write.table(top_1per, file="top_1percent", quote=F ,sep="\t",row.names=F,col.names=F)

# graph the p-values to see whole-genome trends
ggplot(data=top_1per,aes(x=X4, y=X5))+geom_point()+xlab("genome position")+ylab("p-value") + theme_minimal() #+ coord_cartesian(ylim = c(0, max(top_1per$X5)))
# c(0, 200)

# save the plot
ggsave("fishers_plot.pdf")
