#!/usr/bin/env Rscript

#SBATCH --ntasks=4
#SBATCH --mem=60G
#SBATCH --time=10-00:00:00
#SBATCH --job-name='Fst'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/Fst.stdout
#SBATCH -p batch
#SBATCH --array=1
#-8

# how much the file is split up into tasks
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
options(stringsAsFactors = F)

infile <- paste("minor_frequencies",task_id, sep="_")
df <- read_delim(infile,"\t", col_names = FALSE, trim_ws = TRUE)

end <- (task_id*446189‬)
start <- end-446188
outdf <- data.frame()

for (row in start:end){
  p_parent <- df[row,3]
  outdf[row,1:2] <- df[row,1:2]

  for (i in 4:8){
    p_progeny <- df[row,i]
    p_bar <- ((p_parent + p_progeny)/2)
    p_var <- (((p_progeny - p_bar) * (p_progeny - p_bar))/2)
    Fst <- ((p_var)/(p_bar))
    outdf[row,(i-1)] <- Fst
  }
}

colnames(outdf) <- c("chr", "pos", "F18-to-parent", "F27-to-parent", "F28-to-parent", "F50-to-parent", "F58-to-parent")
outfile <- paste("Fst", task_id, sep="_")
write.table(outdf, file=outfile, quote=F ,sep="\t",row.names=F,col.names=T)
