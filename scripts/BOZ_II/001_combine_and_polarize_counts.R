#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --time=1-00:00:00
#SBATCH --job-name='combine_and_polarize'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/BOZ_II/001_combine_and_polarize_counts.stdout
#SBATCH -p koeniglab

# set up the environment
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
options(stringsAsFactors = F)

##### COMBINE PARENT AND PROGENY ALLELE COUNTS
# read in the parent call file to set the chr+pos columns
all_counts <- read_delim("parents.calls","\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
dim(all_counts)

# read in progeny allele count file(s)
# list every file that needs to be read in
file_list <- c("F27A.calls", "F27B.calls", "F50A.calls", "F50B.calls")
F27A <- read_delim(file_list[1],"\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
F27B <- read_delim(file_list[2],"\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
F50A <- read_delim(file_list[3],"\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
F50B <- read_delim(file_list[4],"\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)

# add generation's A & B count sets together
F27 <- F27A[1:2]
F27$ref <- F27A$X3 + F27B$X3
F27$alt <- F27A$X4 + F27B$X4

F50 <- F50A[1:2]
F50$ref <- F50A$X3 + F50B$X3
F50$alt <- F50A$X4 + F50B$X4

# combine parent and progeny counts
all_counts <- cbind(all_counts, F27[,3:4], F50[,3:4])
# add column names for clarity
colnames(all_counts) <- c("CHR","POS","F0_REF", "F0_ALT","F27_REF","F27_ALT","F50_REF","F50_ALT")
# write out the combined count file
write.table(all_counts, file="BOZ_II_counts", quote=F ,sep="\t",row.names=F,col.names=T)


##### POLARIZE PROGENY ALLELES
# calculate polarized progeny AFs if ref allele is minor
# copy the rows with minor ref allele to new df
minor1 <- all_counts[which(all_counts[,3] < all_counts[,4]),]
print(nrow(minor1))
polarized <- c(5, 7)

minor1[,3] <- minor1[, polarized[1]]/(minor1[,5] + minor1[,6])
minor1[,4] <- minor1[, polarized[2]]/(minor1[,7]+ minor1[,8])


# calculate the polarized progeny AFs if alt allele is minor
minor2 <- all_counts[which(all_counts[,4] < all_counts[,3]),]
print(nrow(minor2))
polarized <- c(6, 8)

minor2[,3] <- minor2[, polarized[1]]/(minor2[,5] + minor2[,6])
minor2[,4] <- minor2[, polarized[2]]/(minor2[,7]+ minor2[,8])


# if neither allele is minor (ie. the counts are equal)
equal <- all_counts[which(all_counts[,3] == all_counts[,4]),]
print(nrow(equal))
# I guess just pick a side to be minor?
polarized <- c(5, 7)

equal[,3] <- equal[, polarized[1]]/(equal[,5] + equal[,6])
equal[,4] <- equal[, polarized[2]]/(equal[,7] + equal[,8])


### check that all sites have been accounted for
total <- nrow(minor1)+nrow(minor2)+nrow(equal)
print(paste("total lines:",total))
raw_lines <- nrow(all_counts)
print(paste("line in:",raw_lines))

minor_frequencies <- rbind(minor1[,1:4],minor2[,1:4],equal[,1:4])
# rename columns
colnames(minor_frequencies) <- c("CHR", "POS", "F27_MAF", "F50_MAF")
# add cumulative genome position column
chromosomes <- c("chr1H","chr2H","chr3H","chr4H","chr5H","chr6H","chr7H")
len1 <- sum(minor_frequencies[,1]==chromosomes[1])
len2 <- sum(minor_frequencies[,1]==chromosomes[2]) + len1
len3 <- sum(minor_frequencies[,1]==chromosomes[3]) + len2
len4 <- sum(minor_frequencies[,1]==chromosomes[4]) + len3
len5 <- sum(minor_frequencies[,1]==chromosomes[5]) + len4
len6 <- sum(minor_frequencies[,1]==chromosomes[6]) + len5
len7 <- sum(minor_frequencies[,1]==chromosomes[7]) + len6

for (row in 1:length(minor_frequencies$CHR)) {
  chr_val <- gsub("chr(\\w+)H", "\\1", minor_frequencies[row,1])
    
  if (chr_val=="1"){
    minor_frequencies[row,5] <- minor_frequencies[row,2]
    } else if (chr_val=="2"){
      minor_frequencies[row,5] <- minor_frequencies[row,2] + len1
    } else if (chr_val=="3"){
      minor_frequencies[row,5] <- minor_frequencies[row,2] + len2
    } else if (chr_val=="4"){
      minor_frequencies[row,5] <- minor_frequencies[row,2] + len3
    } else if (chr_val=="5"){
      minor_frequencies[row,5] <- minor_frequencies[row,2] + len4
    } else if (chr_val=="6"){
      minor_frequencies[row,5] <- minor_frequencies[row,2] + len5
    } else if (chr_val=="7"){
      minor_frequencies[row,5] <- minor_frequencies[row,2] + len6
    } else {minor_frequencies[row,5] <- NA
    }
}



write.table(minor_frequencies, file="minor_frequencies_BOZ_II", quote=F ,sep="\t",row.names=F,col.names=T)

