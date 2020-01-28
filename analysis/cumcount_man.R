#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=2-00:00:00
#SBATCH --job-name='cum-man'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/cumulative.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(pacman)
p_load(ggplot2, dplyr, tidyr, data.table)
options(stringsAsFactors = F)
df<-fread("Fisherscumpos")
OutName <- "Fishers"

# parse locus
OutName1<-paste0(OutName, "cumpos")
names(df)[1]<-"CHR"
names(df)[2]<-"POS"
names(df)[17]<-"PVAL"

# format for plotting
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


write.table(result, file=OutName1, quote=F ,sep="\t",row.names=F,col.names=F)

result<-result %>% filter(-log10(result$PVAL)>2)

axisdf = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Manhattan plot
g<-ggplot(result, aes(x=BPcum, y=-log10(PVAL))) +

    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
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
    ylab("-log10(p-value)")

OutName2<-paste0(OutName, "_manhattan.jpeg")
ggsave(OutName2, g, width=10, height=5, units="in")
