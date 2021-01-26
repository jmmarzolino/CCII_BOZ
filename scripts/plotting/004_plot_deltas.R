#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=02:00:00
#SBATCH --job-name='delta_AF'
#SBATCH --output=/rhome/jmarz001/bigdata/CCII_BOZ/scripts/delta_AF.stdout
#SBATCH -p koeniglab


#Set up environment
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")
library(readr)
options(stringsAsFactors = F)

##############################################################################
#####################Read In Dataset##########################################
##############################################################################
#Load in files with change in allele frequencies calculated
delta_AF_F1toALL <- read_delim("delta_AF_F1toALL","\t", col_names = T, trim_ws = TRUE)
delta_AF_DAVIS <- read_delim("delta_AF_DAVIS","\t", col_names = T, trim_ws = TRUE)
delta_AF_BOZ <- read_delim("delta_AF_BOZ","\t", col_names = T, trim_ws = TRUE)

###########################################################################
#Average change in allele frequncy between Davis F18 and Bozeman F27
###########################################################################
#Some plotting to get an idea what the data looks like
##Parents
hist(rowSums(abs(delta_AF_F1toALL[,3])),ylim=c(0,2000000))
hist(rowSums(abs(delta_AF_F1toALL[,4])),ylim=c(0,2000000))
hist(rowSums(abs(delta_AF_F1toALL[,5])),ylim=c(0,2000000))
hist(rowSums(abs(delta_AF_F1toALL[,6])),ylim=c(0,2000000))
hist(rowSums(abs(delta_AF_F1toALL[,7])),ylim=c(0,2000000))
##Davis
hist(rowSums(abs(delta_AF_DAVIS[,3])))
hist(rowSums(abs(delta_AF_DAVIS[,4])))
hist(rowSums(abs(delta_AF_DAVIS[,5])))
##Bozeman
hist(rowSums(abs(delta_AF_BOZ[,3])))
hist(rowSums(abs(delta_AF_BOZ[,4])))
hist(rowSums(abs(delta_AF_BOZ[,5])))

##Mean & median allele frequency change per pop
BOZ_matrix<-data.matrix(delta_AF_BOZ[,3:ncol(delta_AF_BOZ)])
for(i in 1:ncol(BOZ_matrix)){
  print(paste(colnames(BOZ_matrix)[i],"mean is",mean(BOZ_matrix[,i],na.rm=T)))
  print(paste(colnames(BOZ_matrix)[i],"median is",quantile(BOZ_matrix[,i],na.rm=T)[3]))
  print(paste(colnames(BOZ_matrix)[i],"max is",quantile(BOZ_matrix[,i],na.rm=T)[5]))
  print(paste(colnames(BOZ_matrix)[i],"min is",quantile(BOZ_matrix[,i],na.rm=T)[1]))
}
##DAVIS
DAV_matrix<-data.matrix(delta_AF_DAVIS[,3:ncol(delta_AF_DAVIS)])
for(i in 1:ncol(DAV_matrix)){
  print(paste(colnames(DAV_matrix)[i],"mean is",mean(DAV_matrix[,i],na.rm=T)))
  print(paste(colnames(DAV_matrix)[i],"median is",quantile(DAV_matrix[,i],na.rm=T)[3]))
  print(paste(colnames(DAV_matrix)[i],"max is",quantile(DAV_matrix[,i],na.rm=T)[5]))
  print(paste(colnames(DAV_matrix)[i],"min is",quantile(DAV_matrix[,i],na.rm=T)[1]))
}
##PARENTS
F1_matrix<-data.matrix(delta_AF_F1toALL[,3:ncol(delta_AF_F1toALL)])
for(i in 1:ncol(F1_matrix)){
  print(paste(colnames(F1_matrix)[i],"mean is",mean(F1_matrix[,i],na.rm=T)))
  print(paste(colnames(F1_matrix)[i],"median is",quantile(F1_matrix[,i],na.rm=T)[3]))
  print(paste(colnames(F1_matrix)[i],"max is",quantile(F1_matrix[,i],na.rm=T)[5]))
  print(paste(colnames(F1_matrix)[i],"min is",quantile(F1_matrix[,i],na.rm=T)[1]))
}
#could make the print out of that information prettyier/more useful format by putting it into columns or rows and binding those into table/matrix/frame


###########################################################################
#Plot positive and negative delta AF
###########################################################################
#across the genome

#bind all the data columns together
df <- cbind.data.frame(delta_AF_F1toALL,delta_AF_DAVIS[,3:5],delta_AF_BOZ[,3:5])
#write_delim(df,"delta_AF_all",delim="\t",col_names=T)
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
    ylab("Allele Frequency")+
    ylim(0,1)
  
  OutName2<-paste0(OutName, "_manhattan.jpeg")
  ggsave(OutName2, g, width=8, height=5, units="in")
  
  trash <- paste0("X",x)
  names(result)[x]<- trash
}


###########################################################################
#Plot magnitude of delta AF (absolute value) across the genome
###########################################################################
# edit data frame to be absolute values
df <- read_delim("delta_AF_all",delim="\t",col_names=T)
df[,3:ncol(df)] <- abs(df[,3:ncol(df)])

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


###########################################################################
#Plot actual AF changes as distributions
###########################################################################

x=5
gen <- colnames(df)[x]
xlab <- paste("Allele Frequency Change",gen)
OutName <- paste("delta_AF",gen,sep="_")

g <- ggplot(df, aes(get(gen))) + geom_histogram(bins=40)+theme_minimal() + xlab(xlab)+ylab("Frequency")
  
OutName2 <- paste0(OutName, "_distribution.jpeg")
ggsave(g,OutName2, width=10, height=6, units="in")


###########################################################################
#Plot absolute AF changes as distributions
###########################################################################



