#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=2-00:00:00
#SBATCH --job-name='cum-bind'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/sort_combine_threshhold.stdout
#SBATCH -p koeniglab

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(pacman)
p_load(ggplot2, dplyr, tidyr, data.table)
options(stringsAsFactors = F)
pvals <- read_delim("pvals","\t", col_names = FALSE, trim_ws = TRUE)

# Manhattan
# Bonferroni threshold
threshold <- 0.05/nrow(df)
# parse locus
names(df)[1]<-"CHR"

# following code adapted from:
# https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html
# format for plotting
df$BP<-as.numeric(as.character(df[2]))

result <- df %>%
  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(df, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)

result<-result %>% filter(-log10(p_lrt)>2)

axisdf = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Manhattan plot
g<-ggplot(result, aes(x=BPcum, y=-log10(p_lrt))) +

    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
    geom_hline(aes(yintercept=-log10(threshold)), color = "lightslateblue", linetype="dashed", alpha=0.7) +
    scale_color_manual(values = rep(c("black", "grey"), 22 )) +

    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
    scale_y_continuous(expand = c(0, 0.5)) +     # remove space between plot area and x axis

    # Customize the theme:
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
