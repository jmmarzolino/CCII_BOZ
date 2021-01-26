#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=2-00:00:00
#SBATCH --job-name='sort_combine'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/sort_combine_threshhold.stdout
#SBATCH -p batch

# first run the script `cum_pos.R` with a bunch of fixes for multiple file names and different column numbers
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
all_freqs <- read_delim("all_freqs","\t", col_names = FALSE, trim_ws = TRUE)
pvals <- read_delim("pvals","\t", col_names = FALSE, trim_ws = TRUE)
files <- list(all_freqs, pvals)
outnames <- c("freqs_cum_pos", "pval_cum_pos")
chromosomes <- c("chr1H","chr2H","chr3H","chr4H","chr5H","chr6H","chr7H")

for (p in 1:length(files)){
  sample <- data.frame(files[p])
  newcol <- ncol(sample)+1

  len1 <- max(sample[which(sample[,1]==chromosomes[1]),sample[,2]])
  len2 <- max(sample[which(sample[,1]==chromosomes[2]),sample[,2]]) + len1
  len3 <- max(sample[which(sample[,1]==chromosomes[3]),sample[,2]]) + len2
  len4 <- max(sample[which(sample[,1]==chromosomes[4]),sample[,2]]) + len3
  len5 <- max(sample[which(sample[,1]==chromosomes[5]),sample[,2]]) + len4
  len6 <- max(sample[which(sample[,1]==chromosomes[6]),sample[,2]]) + len5
  len7 <- max(sample[which(sample[,1]==chromosomes[7]),sample[,2]]) + len6

  for (row in 1:nrow(sample)) {
    chr_val <- gsub("chr(\\w+)H", "\\1", sample[row,1])

    if (chr_val=="1"){
        sample[row,newcol] <- sample[row,2]
      }
    if (chr_val=="2"){
        sample[row,newcol] <- sample[row,2] + len1
      }
    if (chr_val=="3"){
        sample[row,newcol] <- sample[row,2] + len2
      }
    if (chr_val=="4"){
        sample[row,newcol] <- sample[row,2] + len3
      }
    if (chr_val=="5"){
        sample[row,newcol] <- sample[row,2] + len4
      }
    if (chr_val=="6"){
        sample[row,newcol] <- sample[row,2] + len5
      }
    if (chr_val=="7"){
        sample[row,newcol] <- sample[row,2] + len6
      }}

  # sort the data frames by their cumulative count before writing the table out
  outsample <- sample[order(sample[,newcol]),]
  write.table(outsample, file=outnames[p], quote=F ,sep="\t",row.names=F,col.names=F)
    }

freqs_cum_pos <- read_delim("freqs_cum_pos","\t", col_names = FALSE, trim_ws = TRUE)
pval_cum_pos <- read_delim("pval_cum_pos","\t", col_names = FALSE, trim_ws = TRUE)
# combine the two sorted frames, since they were sorted by the same column, they should be in the same position order
# freqs_cum_pos, pvals
colnames(freqs_cum_pos) <- c("chromosome","position","parents","F18","F27", "F28", "F50", "F58", "cum_position")
colnames(pval_cum_pos) <- c("chromosome","position","pvalue","cum_position")

bind_frame <- cbind(freqs_cum_pos, pval_cum_pos)
write.table(bind_frame, file="bind_frame", quote=F ,sep="\t",row.names=F,col.names=T)

merge_frame <- merge(freqs_cum_pos, pval_cum_pos, by="cum_position")
write.table(merge_frame, file="merge_frame", quote=F ,sep="\t",row.names=F,col.names=T)

# did it work? check if the cumulative positions match by looking for where they don't
bind_mismatch <- bind_frame[which(bind_frame$cum_position != bind_frame$cum_position.1),]
write.table(mismatch, file="cum_pos_bind_mismatch", quote=F ,sep="\t",row.names=F,col.names=T)

merge_mismatch <- merge_frame[which(merge_frame$cum_position != merge_frame$cum_position.1),]
write.table(merge_frame, file="cum_pos_merge_mismatch", quote=F ,sep="\t",row.names=F,col.names=T)

# then parse by the p-value
# set the cutoff for p-values at 1%
cutoff <- quantile(merge_frame$pvalue,0.01)

# if p-value is greater than the cutoff value, replace p-value with NA
# subset and save the data that is at or below the cutoff
merge_frame[which(merge_frame$pvalue >= cutoff),3] <- NA

# convert p-vales into positive values and make 1's into 0's with a -log transform!
merge_frame$transform <- -log10(merge_frame$pvalue)
# save the subsetted table
write.table(merge_frame, file="top_1percent", quote=F ,sep="\t",row.names=F,col.names=T)

# graph the p-values to see whole-genome trends
library(ggplot2)
ggplot(data=merge_frame,aes(x=cum_position, y=transform,color=chromosome.x))+geom_point()+xlab("genome position")+ylab("p-value") + theme_minimal() #+ coord_cartesian(ylim = c(0, max(top_1per$X5)))  # c(0, 200)
# save the plot
ggsave("sort_combine_threshhold_graph.pdf")
