#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --time=01:00:00
#SBATCH --job-name='subset and graph'
#SBATCH --output=plot_top1per_pvals.stdout
#SBATCH -p koeniglab

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(pacman)
p_load(ggplot2, dplyr, tidyr, data.table)

# load table
df <- read_delim("pvals_trim", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
OutName <- "fishers_top1per_binomial"
# set the cutoff for p-values at 1%
cutoff <- quantile(df$X3,0.01,na.rm=T)
# subset and save the data that is at or below the cutoff
top_1per <- df[which(df$X3 <= cutoff),]
# convert p-vales into positive values with a -log transform!
top_1per$X5 <- -log10(top_1per$X3)
# save the subsetted table by writing out
write.table(top_1per, file=OutName, quote=F ,sep="\t",row.names=F,col.names=F)

# set up for plotting with Bonferroni threshold
threshold <- 0.05/nrow(top_1per)
# parse locus
names(top_1per)[1]<-"CHR"
names(top_1per)[2]<-"POS"
top_1per$BP<-as.numeric(top_1per$POS)

result <- top_1per %>%
  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(top_1per, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)

result<-result %>% filter(-log10(top_1per$X3)>2)

axistop_1per = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Manhattan plot
g<-ggplot(result, aes(x=BPcum, y=-log10(top_1per$X3))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
    geom_hline(aes(yintercept=-log10(threshold)), color = "lightslateblue", linetype="dashed", alpha=0.7) +
    scale_color_manual(values = rep(c("black", "grey"), 22 )) +
    # custom X axis:
    scale_x_continuous(label = axistop_1per$CHR, breaks= axistop_1per$center) +
    scale_y_continuous(expand = c(0, 0.5)) +     # remove space between plot area and x axis
    # Custom the theme:
    theme_classic() +
    theme(legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      text=element_text(size=16)) +
    xlab("Chromosome") +
    ylab(expression(-log[10](italic(p))))

OutName1<-paste0(OutName, "_manhattan.jpeg")
ggsave(OutName1, g, width=10, height=5, units="in")

