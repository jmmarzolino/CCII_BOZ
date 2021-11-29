#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --job-name='AF_genome'
#SBATCH --output=AF_genome.stdout
#SBATCH -p short

options(stringsAsFactors = F)
library(pacman)
library(readr)
p_load(ggplot2, dplyr, tidyr, data.table)

# set variables
FileName="minor_frequencies"

# import data
df <- read_delim(FileName,"\t", col_names = FALSE, trim_ws = TRUE)

# parse locus
names(df)[1]<-"CHR"
names(df)[2]<-"POS"
col_names <- c("CHR","POS","F1","F18", "F27","F28","F50","F58")
#names(df)[10]<-"BPcum"

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

axisresult = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

for (x in 3:8){
  OutName <- paste0("AF_genome",col_names[x])
  names(result)[x]<-"Y"
  xlab <- col_names[x]

  g <- ggplot(result, aes(x=BPcum, y=Y)) +
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
  scale_color_manual(values = rep(c("black", "grey"), 22 )) +
  # custom X axis:
  scale_x_continuous(label = axisresult$CHR, breaks= axisresult$center) +
  scale_y_continuous(expand = c(0, 0.5)) +     # remove space between plot area and x axis
  # Customize the theme:
  theme_classic() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text=element_text(size=16)) +
  xlab(xlab) +
  ylab("Allele Frequency")+
  ylim(0,1)

  OutName2<-paste0(OutName, "_manhattan.jpeg")
  ggsave(OutName2, g, width=8, height=5, units="in")

  trash <- paste0("X",x)
  names(result)[x]<- trash
}

