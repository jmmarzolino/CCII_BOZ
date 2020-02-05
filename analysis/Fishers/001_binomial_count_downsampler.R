library(readr)
library(stats)
df <- read_delim("filtered_counts", "\t", col_names = FALSE,trim_ws = TRUE)

binom_counts <- df[1:2]
rows <- nrow(df)

Dav_mean_counts <- round(((df$X7+df$X8+df$X11+df$X12+df$X15+df$X16)/3),digits=0)

Boz27_ref_prob <- (df$X9)/(df$X9+df$X10)
Boz50_ref_prob <- (df$X13)/(df$X13+df$X14)
parent_ref_prob <- (df$X5)/(df$X5+df$X6)

Boz27_binom_ref_count <- rbinom(n=rows,size=Dav_mean_counts,prob=Boz27_ref_prob)
parent_binom_ref_count <- rbinom(n=rows,size=Dav_mean_counts,prob=parent_ref_prob)
Boz50_binom_ref_count <- rbinom(n=rows,size=Dav_mean_counts,prob=Boz50_ref_prob)

# add the new parental ref counts to the out df
binom_counts[,3] <- parent_binom_ref_count
# use the ref counts and Dav mean counts to create the alt  counts
binom_counts[,4] <- (Dav_mean_counts)-(parent_binom_ref_count)
# copy over Davis counts F18
binom_counts[,5:6] <- df[,7:8]

# repeat above with Boz populations
binom_counts[,7] <- Boz27_binom_ref_count
binom_counts[,8] <- (Dav_mean_counts)-(Boz27_binom_ref_count)
# copy Davis F28
binom_counts[,9:10] <- df[,11:12]

# Boz F50
binom_counts[,11] <- Boz50_binom_ref_count
binom_counts[,12] <- (Dav_mean_counts)-(Boz50_binom_ref_count)
# Davis F58
binom_counts[,13:14] <- df[,15:16]

write_delim(binom_counts,"binomial_counts",delim="\t",col_names = F)
