#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --job-name='man plot'
#SBATCH --output=man_plot.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(pacman)
p_load(ggplot2, dplyr, tidyr, data.table)
options(stringsAsFactors = F)

df<-fread("Fst_cumcumpos")

# parse locus
names(df)[1]<-"CHR"
names(df)[2]<-"POS"
col_names <- c("CHR","POS","F18-to-parent", "F27-to-parent", "F28-to-parent", "F50-to-parent", "F58-to-parent")
names(df)[10]<-"BPcum"

for (x in 3:7){
  OutName <- paste0("Fst",x)
  names(df)[x]<-"Y"
  xlab <- col_names[x]
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
  ylim(NA,0.25)

  OutName2<-paste0(OutName, "_manhattan.jpeg")
  ggsave(OutName2, g, width=8, height=5, units="in")

  trash <- paste0("X",x)
  names(df)[x]<- trash
}
