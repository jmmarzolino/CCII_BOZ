#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=01:00:00
#SBATCH --job-name='Fst_genome_plot'
#SBATCH --output=/bigdata/koeniglab/jmarz001/CCII_BOZ/scripts/Fst_genome_plot.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(pacman)
p_load(ggplot2, dplyr, tidyr, data.table)
options(stringsAsFactors = F)

df<-read_delim("Fst_cumpos","\t", col_names = T, trim_ws = TRUE)

for (x in 3:7){
  # record the columns name for graph
  col <- colnames(df)[x]
  xlab <- col
  OutName <- paste("Fst",col,sep="_")
  # update working column name to Y for ggplotting
  names(df)[x]<-"Y"

  axisdf = df %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

  g <- ggplot(df, aes(x=BPcum, y=Y)) +
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
  xlab(xlab) +
  ylab("Fst")+
  ylim(0,1)

  OutName2<-paste0(OutName, "_manhattan.jpeg")
  ggsave(OutName2, g, width=8, height=5, units="in")

  names(df)[x]<- paste0("X",x)
}

