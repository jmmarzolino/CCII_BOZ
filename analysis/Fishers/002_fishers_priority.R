#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --job-name='fishers'
#SBATCH --output=fishers_priority.stdout
#SBATCH -p koeniglab
#SBATCH --time=5-00:00:00
#SBATCH --array=1-13

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
progeny_counts_filtered <- read_delim("binomial_counts",
    "\t", escape_double = FALSE, col_names = FALSE,
    trim_ws = TRUE)
df <- na.omit(progeny_counts_filtered)

end <- (274429*task_id‬)
start <- (end-274428)
# now compute the p-values
for (row in start:end) {
  # take allele counts from each line in turn, convert into a matrix with generations as rows
  site_counts <- matrix(as.numeric(df[row,3:14]),nrow=6,byrow=T)
  if (sum(site_counts)!=0){
    p_value <- fisher.test(site_counts,workspace = 2e8, hybrid=T)[1]
    df[row,15] <- p_value
  }
}
file_name <- paste("filtered_Fishers_pvalues_allpops", sep="", task_id)
write.table(df[start:end,], file=file_name, quote=F ,sep="\t",row.names=F,col.names=F)
