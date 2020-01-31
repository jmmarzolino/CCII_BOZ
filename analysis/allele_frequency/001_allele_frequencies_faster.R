#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --job-name='AF'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/delta_AF.stdout
#SBATCH -p batch

# set up environment
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
options(stringsAsFactors = F)

# read in allele count file
progeny_counts <- read_delim("filtered_counts","\t", col_names = FALSE, trim_ws = TRUE)

# calculate polarized progeny AFs if ref allele is minor
# copy the rows with minor ref allele to new df
minor1 <- progeny_counts[which(progeny_counts[,5] < progeny_counts[,6]),]
print(nrow(minor1))
polarized <- c(7,9,11,13,15)

minor1[,3] <- minor1[,5]/(minor1[,5] + minor1[,6])
minor1[,4] <- minor1[,polarized[1]]/(minor1[,7]+ minor1[,8])
minor1[,5] <- minor1[,polarized[2]]/(minor1[,9]+ minor1[,10])
minor1[,6] <- minor1[,polarized[3]]/(minor1[,11]+ minor1[,12])
minor1[,7] <- minor1[,polarized[4]]/(minor1[,13]+ minor1[,14])
minor1[,8] <- minor1[,polarized[5]]/(minor1[,15]+ minor1[,16])


# calculate the polarized progeny AFs if alt allele is minor
minor2 <- progeny_counts[which(progeny_counts[,6] < progeny_counts[,5]),]
print(nrow(minor2))
polarized <- c(8,10,12,14,16)

minor2[,3] <- minor2[,6]/(minor2[,5] + minor2[,6])
minor2[,4] <- minor2[,polarized[1]]/(minor2[,7]+ minor2[,8])
minor2[,5] <- minor2[,polarized[2]]/(minor2[,9]+ minor2[,10])
minor2[,6] <- minor2[,polarized[3]]/(minor2[,11]+ minor2[,12])
minor2[,7] <- minor2[,polarized[4]]/(minor2[,13]+ minor2[,14])
minor2[,8] <- minor2[,polarized[5]]/(minor2[,15]+ minor2[,16])

# if neither allele is minor (ie. the counts are equal)
equal <- progeny_counts[which(progeny_counts[,5] == progeny_counts[,6]),]
print(nrow(equal))
# I guess just pick a side to be minor?
polarized <- c(7,9,11,13,15)

equal[,3] <- equal[,5]/(equal[,5] + equal[,6])
equal[,4] <- equal[,polarized[1]]/(equal[,7]+ equal[,8])
equal[,5] <- equal[,polarized[2]]/(equal[,9]+ equal[,10])
equal[,6] <- equal[,polarized[3]]/(equal[,11]+ equal[,12])
equal[,7] <- equal[,polarized[4]]/(equal[,13]+ equal[,14])
equal[,8] <- equal[,polarized[5]]/(equal[,15]+ equal[,16])

total <- nrow(minor1)+nrow(minor2)+nrow(equal)
print(paste("total lines:",total))
raw_lines <- nrow(progeny_counts)
print(paste("line in:",raw_lines))

minor_frequencies <- rbind(minor1[,1:8],minor2[,1:8],equal[,1:8])
write.table(minor_frequencies, file="minor_frequencies", quote=F ,sep="\t",row.names=F,col.names=F)
