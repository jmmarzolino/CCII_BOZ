#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --job-name='fishers'
#SBATCH -p koeniglab
#SBATCH --time=5-00:00:00
#SBATCH --array=1-1690‬
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/fishers_priority_{$SLURM_ARRAY_TASK_ID}.stdout

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
progeny_counts_filtered <- read_delim("filtered_counts",
    "\t", escape_double = FALSE, col_names = FALSE,
    trim_ws = TRUE)

# define rows in the matrix to grab for pairwise (and bigger) comparisons
parents <- 1
F18 <- 2
F27 <- 3
F28 <- 4
F50 <- 5
F58 <- 6

#compare everything
all6pops <- c(parents, F18, F27, F28, F50, F58)

end <- (2111*task_id‬)
start <- (end-2110)

# now compute the p-values
for (row in start:end) {
  # take allele counts from each line in turn, convert into a matrix with generations as rows
  site_counts <- matrix(as.numeric(progeny_counts_filtered[row,5:16]),nrow=6,byrow=T)

  p_value <- fisher.test(site_counts[all6pops,1:2],workspace = 2e8, hybrid=T)[1]
  progeny_counts_filtered[row,17] <- p_value
}
file_name <- paste("filtered_Fishers_pvalues_allpops", sep="", task_id)
write.table(progeny_counts_filtered[start:end,], file=file_name, quote=F ,sep="\t",row.names=F,col.names=F)
