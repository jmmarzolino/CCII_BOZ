#!/usr/bin/env Rscript

#SBATCH --ntasks=4
#SBATCH --mem=60G
#SBATCH --time=10-00:00:00
#SBATCH --job-name='AF'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/delta_AF.stdout
#SBATCH -p koeniglab
#SBATCH --array=1-8

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
options(stringsAsFactors = F)

progeny_counts <- read_delim("progeny_counts","\t", col_names = FALSE, trim_ws = TRUE)

# how much the file is split up into tasks
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

end <- (task_id*446189‬)
start <- end-446188
minor_frequencies <- data.frame()

for (row in start:end){
  minor_frequencies[row,1] <- progeny_counts[row,1]
  minor_frequencies[row,2] <- progeny_counts[row,2]

  x1 <- progeny_counts[row,5]
  x2 <- progeny_counts[row,6]

  if (x1 < x2){
    minor <- x1
    polarized <- c(7,9,11,13,15)
  }
  else if (x2 < x1) {
    minor <- x2
    #base <- progeny_counts[row,4]
    polarized <- c(8,10,12,14,16)
  }
  # calculate the minor & polarized progeny AFs
  minor_frequencies[row,3] <- minor/(x1 + x2)
  minor_frequencies[row,4] <- progeny_counts[row,polarized[1]]/(progeny_counts[row,7]+ progeny_counts[row,8])
  minor_frequencies[row,5] <- progeny_counts[row,polarized[2]]/(progeny_counts[row,9]+ progeny_counts[row,10])
  minor_frequencies[row,6] <- progeny_counts[row,polarized[3]]/(progeny_counts[row,11]+ progeny_counts[row,12])
  minor_frequencies[row,7] <- progeny_counts[row,polarized[4]]/(progeny_counts[row,13]+ progeny_counts[row,14])
  minor_frequencies[row,8] <- progeny_counts[row,polarized[5]]/(progeny_counts[row,15]+ progeny_counts[row,16])
}
outfilename <- paste("minor_frequencies",task_id, sep="_")

write.table(minor_frequencies, file=outfilename, quote=F ,sep="\t",row.names=F,col.names=F)
