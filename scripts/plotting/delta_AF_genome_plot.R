#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=02:00:00
#SBATCH --job-name='delta_AF'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/plot_delta_AF_genome.stdout
#SBATCH -p koeniglab

#Set up environment
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
library(pacman)
p_load(ggplot2, dplyr, tidyr, data.table)
options(stringsAsFactors = F)

##############################################################################
#####################Read In Dataset##########################################
##############################################################################
#Load in files with change in allele frequencies calculated
delta_AF_F1toALL <- read_delim("delta_AF_F1toALL","\t", col_names = T, trim_ws = TRUE)
delta_AF_DAVIS <- read_delim("delta_AF_DAVIS","\t", col_names = T, trim_ws = TRUE)
delta_AF_BOZ <- read_delim("delta_AF_BOZ","\t", col_names = T, trim_ws = TRUE)


###########################################################################
#Plot positive and negative delta AF across the genome
###########################################################################

#bind all the data columns together
df <- cbind.data.frame(delta_AF_F1toALL,delta_AF_DAVIS[,3:5],delta_AF_BOZ[,3:5])
# name relevant columns
names(df)[1]<-"CHR"
names(df)[2]<-"POS"

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

#head(result)
axisresult = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

for (x in 3:13){
  OutName <- paste0(colnames(result)[x],"_deltaAF")
  xlab <- colnames(result)[x]
  names(result)[x]<-"Y"

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
    ylab("Change in Allele Frequency")+
    ylim(-1,1)

  OutName2<-paste0(OutName, "_manhattan.jpeg")
  ggsave(OutName2, g, width=10, height=6, units="in")

  trash <- paste0("X",x)
  names(result)[x]<- trash
}


###########################################################################
#Plot magnitude of delta AF (absolute value) across the genome
###########################################################################
# edit data frame to be absolute values
df <- cbind.data.frame(abs(delta_AF_F1toALL,delta_AF_DAVIS[,3:5],delta_AF_BOZ[,3:5]))

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

#head(result)
axisresult = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

for (x in 3:13){
  OutName <- paste0(colnames(result)[x],"_deltaAF_abs")
  xlab <- colnames(result)[x]
  names(result)[x]<-"Y"

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
    ylab("Change in Allele Frequency")+
    ylim(0,1)

  OutName2<-paste0(OutName, "_manhattan.jpeg")
  ggsave(OutName2, g, width=10, height=6, units="in")

  trash <- paste0("X",x)
  names(result)[x]<- trash
}

