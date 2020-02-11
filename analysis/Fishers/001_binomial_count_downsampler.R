#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=10:00:00
#SBATCH --job-name='binomial_counts'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/binomial_count_downsampler.stdout
#SBATCH -p koeniglab

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
options(stringsAsFactors = F)
library(readr)
library(stats)
df <- read_delim("filtered_counts", "\t", col_names = FALSE,trim_ws = TRUE)

# copy over chr & pos columns, parental & F18 Davis counts
binom_counts <- df[,1:2] #initialize binom_counts df first
binom_counts[,3:6] <- df[,5:8]
rows <- nrow(df)

Dav_mean_counts <- (df$X7+df$X8+df$X11+df$X12+df$X15+df$X16)/3
Dav_mean_counts <- round(Dav_mean_counts,digits=0)

Boz27_ref_prob <- (df$X9)/(df$X9+df$X10)
Boz50_ref_prob <- (df$X13)/(df$X13+df$X14)

Boz27_binom_ref_count <- rbinom(n=rows,size=Dav_mean_counts,prob=Boz27_ref_prob)
Boz50_binom_ref_count <- rbinom(n=rows,size=Dav_mean_counts,prob=Boz50_ref_prob)

# add the new BOZ ref counts to the out df
binom_counts[,7] <- Boz27_binom_ref_count
# use the ref counts and Dav mean counts to create the alt  counts
binom_counts[,8] <- (Dav_mean_counts)-(Boz27_binom_ref_count)
# copy Davis F28
binom_counts[,9:10] <- df[,11:12]

# Boz F50
binom_counts[,11] <- Boz50_binom_ref_count
binom_counts[,12] <- (Dav_mean_counts)-(Boz50_binom_ref_count)
# Davis F58
binom_counts[,13:14] <- df[,15:16]

write_delim(binom_counts,"binomial_counts",delim="\t",col_names = F)
