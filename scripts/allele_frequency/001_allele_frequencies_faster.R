#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --job-name='AF'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/AF002_allele_frequencies.stdout
#SBATCH -p batch

# set up environment
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
options(stringsAsFactors = F)

# read in allele count file
counts <- read_delim("filtered_counts","\t", col_names = T, trim_ws = TRUE)

# calculate polarized progeny AFs if ref allele is minor
# copy the rows with minor ref allele to new df
REF_minor <- counts[which(counts[,3] < counts[,4]),]
print(nrow(REF_minor))
polarized <- c(5,7,9,11,13)
minor1 <- REF_minor[,1:2]

minor1[,3] <- REF_minor[,3]/(REF_minor[,3] + REF_minor[,4])
minor1[,4] <- REF_minor[,polarized[1]]/(REF_minor[,5]+ REF_minor[,6])
minor1[,5] <- REF_minor[,polarized[2]]/(REF_minor[,7]+ REF_minor[,8])
minor1[,6] <- REF_minor[,polarized[3]]/(REF_minor[,9]+ REF_minor[,10])
minor1[,7] <- REF_minor[,polarized[4]]/(REF_minor[,11]+ REF_minor[,12])
minor1[,8] <- REF_minor[,polarized[5]]/(REF_minor[,13]+ REF_minor[,14])
colnames(minor1) <- c("CHR","POS","F0 MAF","F18 MAF","F27 MAF","F28 MAF","F50 MAF","F58 MAF")

# calculate the polarized progeny AFs if alt allele is minor
ALT_minor <- counts[which(counts[,4] < counts[,3]),]
print(nrow(ALT_minor))
polarized <- c(6,8,10,12,14)
minor2 <- ALT_minor[,1:2]

minor2[,3] <- ALT_minor[,4]/(ALT_minor[,3] + ALT_minor[,4])
minor2[,4] <- ALT_minor[,polarized[1]]/(ALT_minor[,5]+ ALT_minor[,6])
minor2[,5] <- ALT_minor[,polarized[2]]/(ALT_minor[,7]+ ALT_minor[,8])
minor2[,6] <- ALT_minor[,polarized[3]]/(ALT_minor[,9]+ ALT_minor[,10])
minor2[,7] <- ALT_minor[,polarized[4]]/(ALT_minor[,11]+ ALT_minor[,12])
minor2[,8] <- ALT_minor[,polarized[5]]/(ALT_minor[,13]+ ALT_minor[,14])
colnames(minor2) <- c("CHR","POS","F0 MAF","F18 MAF","F27 MAF","F28 MAF","F50 MAF","F58 MAF")

# if neither allele is minor (ie. the counts are equal)
equal <- counts[which(counts[,3] == counts[,4]),]
print(nrow(equal))
# I guess just pick a side to be minor?
polarized <- c(5,7,9,11,13)
minor3 <- equal[,1:2]

minor3[,3] <- equal[,3]/(equal[,3] + equal[,4])
minor3[,4] <- equal[,polarized[1]]/(equal[,5]+ equal[,6])
minor3[,5] <- equal[,polarized[2]]/(equal[,7]+ equal[,8])
minor3[,6] <- equal[,polarized[3]]/(equal[,9]+ equal[,10])
minor3[,7] <- equal[,polarized[4]]/(equal[,11]+ equal[,12])
minor3[,8] <- equal[,polarized[5]]/(equal[,13]+ equal[,14])
colnames(minor3) <- c("CHR","POS","F0 MAF","F18 MAF","F27 MAF","F28 MAF","F50 MAF","F58 MAF")

total <- nrow(minor1)+nrow(minor2)+nrow(minor3)
print(paste("total lines:",total))
raw_lines <- nrow(counts)
print(paste("lines in:",raw_lines))

minor_frequencies <- rbind(minor1,minor2,minor3)
write.table(minor_frequencies, file="minor_frequencies", quote=F ,sep="\t",row.names=F,col.names=T)
