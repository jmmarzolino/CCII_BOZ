#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --job-name='delta He'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/delta_He.stdout
#SBATCH -p koeniglab

# set up environment
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(ggplot2)
options(stringsAsFactors = F)

# list the generations to name files
gens <- c("F0","F18","F27","F28","F50","F58")

i=1
# read in parental He file to use as all_pop file base
OutName <- paste(gens[i],"expected_het",sep="_")
all_He <- read_delim(OutName, delim="\t", col_names = T, trim_ws = T)
colnames(all_He)[i+2] <- gens[i]

all_weighted_He <- all_He[,c(1,2,4)]
colnames(all_weighted_He)[i+2] <- gens[i]

for(i in 2:length(gens)){
  # read in each pop file
  OutName <- paste(gens[i],"expected_het",sep="_")
  df <- read_delim(OutName, delim="\t", col_names = T, trim_ws = T)
  # copy He column from progeny as new, pop-named He column
  all_He[i+2] <- df[,3]
  colnames(all_He)[i+2] <- gens[i]
  # copy weighted He column too
  all_weighted_He[i+2] <- df[,4]
  colnames(all_weighted_He)[i+2] <- gens[i]
}

# write out the all population file!
write_delim(all_He,"all_He",delim="\t",col_names = T)
write_delim(all_weighted_He,"all_weighted_He",delim="\t",col_names = T)
################################################################################
######################  DELTA   ################################################
################################################################################
# initialize a data frame to write to
delta_He <- all_He[1:2]
# calculate difference in He between all progeny and parents
delta_He$F0toF18 <- all_He$F18 - all_He$F0
delta_He$F0toF28 <- all_He$F28 - all_He$F0
delta_He$F0toF58 <- all_He$F58 - all_He$F0
delta_He$F0toF27 <- all_He$F27 - all_He$F0
delta_He$F0toF50 <- all_He$F50 - all_He$F0
write_delim(delta_He,"delta_He_F0toALL",delim="\t",col_names = T)

# between each Davis generation
delta_He <- all_He[,1:2]
delta_He$F18toF28 <- all_He$F28 - all_He$F18
delta_He$F28toF58 <- all_He$F58 - all_He$F28
delta_He$F18toF58 <- all_He$F58 - all_He$F18
write_delim(delta_He,"delta_He_DAVIS",delim="\t",col_names = T)

# between each Bozeman generation
delta_He <- all_He[,1:2]
delta_He$F18toF27 <- all_He$F27 - all_He$F18
delta_He$F27toF50 <- all_He$F50 - all_He$F27
delta_He$F18toF50 <- all_He$F50 - all_He$F18
write_delim(delta_He,"delta_He_BOZ",delim="\t",col_names = T)

######################  WEIGHTED  ##############################################
################################################################################
# initialize a data frame to write to
delta_weighted_He <- all_weighted_He[1:2]
# calculate difference in He between all progeny and parents
delta_weighted_He$F0toF18 <- all_weighted_He$F18 - all_weighted_He$F0
delta_weighted_He$F0toF28 <- all_weighted_He$F28 - all_weighted_He$F0
delta_weighted_He$F0toF58 <- all_weighted_He$F58 - all_weighted_He$F0
delta_weighted_He$F0toF27 <- all_weighted_He$F27 - all_weighted_He$F0
delta_weighted_He$F0toF50 <- all_weighted_He$F50 - all_weighted_He$F0
write_delim(delta_weighted_He,"delta_weighted_He_F0toALL",delim="\t",col_names = T)

# between each Davis generation
delta_weighted_He <- all_weighted_He[,1:2]
delta_weighted_He$F18toF28 <- all_weighted_He$F28 - all_weighted_He$F18
delta_weighted_He$F28toF58 <- all_weighted_He$F58 - all_weighted_He$F28
delta_weighted_He$F18toF58 <- all_weighted_He$F58 - all_weighted_He$F18
write_delim(delta_weighted_He,"delta_weighted_He_DAVIS",delim="\t",col_names = T)

# between each Bozeman generation
delta_weighted_He <- all_weighted_He[,1:2]
delta_weighted_He$F18toF27 <- all_weighted_He$F27 - all_weighted_He$F18
delta_weighted_He$F27toF50 <- all_weighted_He$F50 - all_weighted_He$F27
delta_weighted_He$F18toF50 <- all_weighted_He$F50 - all_weighted_He$F18
write_delim(delta_weighted_He,"delta_weighted_He_BOZ",delim="\t",col_names = T)


################################################################################
######################  EXPLORE THE DATA  ######################################
################################################################################

# note that to get full results you need to go back and regenerate each of the data sets that got put into files or load in the files {fix laaaater}
# convert numeric columns to matrix for ease of use
m <- as.matrix(delta_He[3:ncol(delta_He)])
# what's the mean/median/min/max of the delta He and delta weighted He data sets FOR EACH COLUMN
for(i in 1:ncol(m)){
  print(round(mean(m[,i],na.rm = T),digits=3))
}

library(matrixStats)
round(colQuantiles(m,na.rm = T),digits=3)


# convert numeric columns to matrix for ease of use
m <- as.matrix(delta_weighted_He[3:ncol(delta_weighted_He)])
for(i in 1:ncol(m)){
  print(round(mean(m[,i],na.rm = T),digits=3))
}

library(matrixStats)
round(colQuantiles(m,na.rm = T),digits=3)
