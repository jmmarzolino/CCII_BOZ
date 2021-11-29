#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --time=01:00:00
#SBATCH --job-name='fishers_plot'
#SBATCH --output=plot_all_pvals.stdout
#SBATCH -p koeniglab

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(pacman)
p_load(ggplot2, dplyr, tidyr, data.table)

# load table
df <- read_delim("pvals_trim", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
OutName <- "fishers_binomial"

# convert p-vales into positive values with a -log transform!
df$transform <- -log10(df$X3)
# set up for plotting with Bonferroni threshold
threshold <- -log10(0.05/nrow(df))
# parse locus
names(df)[1]<-"CHR"
names(df)[2]<-"POS"
df$BP<-as.numeric(df$POS)

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

result<-result %>% filter(-log10(df$X3)>2)
axisdf = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Manhattan plot
g<-ggplot(result, aes(x=BPcum, y=transform)) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
    geom_hline(aes(yintercept=threshold), color = "lightslateblue", linetype="dashed", alpha=0.7) +
    scale_color_manual(values = rep(c("black", "grey"), 22 )) +
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
    scale_y_continuous(expand = c(0, 0.5)) +     # remove space between plot area and x axis
    # Custom the theme:
    theme_classic() +
    theme(legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      text=element_text(size=18)) +
    xlab("Chromosome") +
    ylab(expression(-log[10](italic(p))))

OutName1<-paste0(OutName, "_genome.png")
ggsave(OutName1, g, width=10, height=5, units="in")
