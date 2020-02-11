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
# read in allele count file
sample <- read_delim("filtered_counts", delim="\t", col_names = FALSE, trim_ws = TRUE)

# list the generations to name out files
gens <- c("F0","F18","F27","F28","F50","F58")
OutName <- paste(gens[task_id],"expected_het",sep="_")

# make a new df for the population's heterozygosity and weight value
df <- sample[1:2]
# assign the pop's column counts to ref/alt variables
ref_count <- sample[,((2*(task_id + 1))+1)]
alt_count <- sample[,(2*(task_id+2))]
# record total number of reads divided by max no. reads in that pop
count_sum <- (ref_count +alt_count)
pop_max <- max(count_sum)
weight <- (count_sum)/(pop_max)
# calculate He
het <- 2*(ref_count/(ref_count + alt_count))*(alt_count/(ref_count + alt_count))
# Weight the expected heterozygosity fraction by count
weighted_het <-het*weight

# add het and weighted het to df
df <- cbind(df, het, weighted_het)
# name the columsn for clarity
colnames(df) <- c("CHR","POS","expected_heterozygosity","weighted_ex_heterozygosity")
# write out the population files with He and the weighted He
write.table(df, file=OutName, quote=F ,sep="\t",row.names=F,col.names=T)
# follow up by plotting for comparison

#df <- read_delim(OutName,sep="\t",col.names=T,trim_ws=T)
# overall contribution of Heterozygosity is weighted by our /certainty/ about genotpye sequencing
binlist <- c(5,10,15, 20, 30, 40, 50)

xlab <- paste(gens[task_id], "Expected Heterozygosity", sep=" ")
for (m in binlist) {
  OutName2 = paste(paste(OutName,m,sep="_"),"bins.pdf",sep="")
  ggplot(df, aes(het)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2750000)
  ggsave(OutName2, width = 10, height = 8, units = "in")
}

xlab <- paste(gens[task_id], "Weighted Expected Heterozygosity", sep=" ")
for (m in binlist) {
  OutName2 = paste(paste(paste(OutName,"weighted",sep="_"),m,sep="_"),"bins.pdf",sep="")
  ggplot(df, aes(weighted_het)) + geom_histogram(bins=m)+ theme_minimal()+xlab(xlab)+ylim(0,2750000)
  ggsave(outname, width = 10, height = 8, units = "in")
}
