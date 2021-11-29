#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=02:00:00
#SBATCH --job-name='delta_AF'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/delta_AF.stdout
#SBATCH -p koeniglab

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
options(stringsAsFactors = F)

minor_frequencies <- read_delim("minor_frequencies","\t", col_names = FALSE, trim_ws = TRUE)
colnames(minor_frequencies) <- c("CHR","POS","F1","F18","F27","F28","F50","F58")
delta_AF <- minor_frequencies[,1:2]

# between all progeny and parents
delta_AF$F1toF18 <- minor_frequencies$F18 - minor_frequencies$F1
delta_AF$F1toF28 <- minor_frequencies$F28 - minor_frequencies$F1
delta_AF$F1toF58 <- minor_frequencies$F58 - minor_frequencies$F1
delta_AF$F1toF27 <- minor_frequencies$F27 - minor_frequencies$F1
delta_AF$F1toF50 <- minor_frequencies$F50 - minor_frequencies$F1
write_delim(delta_AF,"delta_AF_F1toALL",delim="\t",col_names = T)

# between each Davis generation
delta_AF <- minor_frequencies[,1:2]
delta_AF$F18toF28 <- minor_frequencies$F28 - minor_frequencies$F18
delta_AF$F28toF58 <- minor_frequencies$F58 - minor_frequencies$F28
delta_AF$F18toF58 <- minor_frequencies$F58 - minor_frequencies$F18
write_delim(delta_AF,"delta_AF_DAVIS",delim="\t",col_names = T)

# between each Bozeman generation
delta_AF <- minor_frequencies[,1:2]
delta_AF$F18toF27 <- minor_frequencies$F27 - minor_frequencies$F18
delta_AF$F27toF50 <- minor_frequencies$F50 - minor_frequencies$F27
delta_AF$F18toF50 <- minor_frequencies$F50 - minor_frequencies$F18
write_delim(delta_AF,"delta_AF_BOZ",delim="\t",col_names = T)

