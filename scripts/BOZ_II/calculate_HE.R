library(readr)
library(ggplot2)
setwd("/rhome/jmarz001/bigdata/CCII_BOZ/results")

F27A <- read_delim("F27A.calls", "\t", col_names = F)
F27B <- read_delim("F27B.calls", "\t", col_names = F)
F50A <- read_delim("F50A.calls", "\t", col_names = F)
F50B <- read_delim("F50A.calls", "\t", col_names = F)

F27 <- F27A[,1:2]
F27$X3 <- F27A$X3 + F27B$X3
F27$X4 <- F27A$X4 + F27B$X4


F50 <- F50A[,1:2]
F50$X3 <- F50A$X3 + F50B$X3
F50$X4 <- F50A$X4 + F50B$X4

### ~very loosely~ calculate coverage across the parental sites you have ##
(sum(F27$X3) + sum(F27$X4))/nrow(F27)
(sum(F50$X3) + sum(F50$X4))/nrow(F50)
# coverage for F50 is a bit higher than for F27

#### CALCULATE ALLELE FREQUENCIES #######
F50$p <- F50$X3 / (F50$X3+F50$X4)
F50$q <- 1 - F50$p

F27$p <- F27$X3 / (F27$X3+F27$X4)
F27$q <- 1 - F27$p

######## CALCULATE EXPECTED HETEROZYGOSITY #############
F50$heterozygosity <- 2* (F50$X3/ (F50$X3+ F50$X4)) *  (F50$X4/ (F50$X3+ F50$X4))
F27$heterozygosity <- 2* (F27$X3/ (F27$X3+ F27$X4)) *  (F27$X4/ (F27$X3+ F27$X4))

### WRITE OUT CALCULATIONS ###
write.table(F50, file="F50_het.calls", quote=F ,sep="\t",row.names=F,col.names=F)
write.table(F27, file="F27_het.calls", quote=F ,sep="\t",row.names=F,col.names=F)


##### PLOT THE HETEROZYGOSITY RELATIONSHIP ####
mean(F50$heterozygosity, na.rm=T)
mean(F27$heterozygosity, na.rm=T)

generation <- c(27, 50)
heterozygosity <- c(mean(F27$heterozygosity, na.rm=T), mean(F50$heterozygosity, na.rm=T))

het_decay <- data.frame(generation, heterozygosity)

ggplot(het_decay, aes(x=generation, y=heterozygosity)) + 
  geom_line(aes()) + 
  geom_point(aes()) + 
  theme_bw()



## He for all generations from previous sequencing, for reference
#parents	0.2048737
#F18 	0.1718248
#F27 	0.108593
#F28 	0.147064
#F50 	0.1740904
#F58 	0.058337

CCII <- c("Davis", "Davis", "Davis", "Davis", "Bozeman", "Bozeman")
generation <- c(0, 18, 28, 58, 27, 50)
heterozygosity1 <- c(0.2048737, 0.1718248, 0.147064, 0.058337, 0.108593, 0.1740904)
heterozygosity2 <- c((0.2048737/2), (0.1718248/2), (0.147064/2), (0.058337/2), mean(F27$heterozygosity, na.rm=T), mean(F50$heterozygosity, na.rm=T))

het_decay <- data.frame(CCII,generation, heterozygosity1, heterozygosity2)

ggplot(het_decay, aes(x=generation, y=heterozygosity1, group=CCII)) +
  geom_line(aes(linetype=CCII)) +
  geom_point(aes(shape=CCII)) +
  theme_bw()

ggplot(het_decay, aes(x=generation, y=heterozygosity2, group=CCII)) +
  geom_line(aes(linetype=CCII)) +
  geom_point(aes(shape=CCII)) +
  theme_bw()


### compare number of covered sites ###
sum((F27$X3+F27$X4)>0)
sum((F50$X3+F50$X4)>0)
nrow(F27)
sum((F27$X3+F27$X4)>0) / nrow(F27)
sum((F50$X3+F50$X4)>0) / nrow(F50)


### compare number of segregating sites ###
sum(F27$X3 != 0 & F27$X4 !=0)/nrow(F27)
sum(F50$X3 != 0 & F50$X4 !=0)/nrow(F50)


######################### ARE COUNTS/COVERAGES SIG. DIFFERENT? ##############################
merge <- cbind(F27[,1:4], F50[,3:4])
fsh_mtrx <- matrix(as.numeric(merge[1,3:6]), nrow=2,byrow = TRUE)
fisher.test(fsh_mtrx)


## Allele Frequencies
F27ref <- (sum(F27$X3)/(sum(F27$X3) + sum(F27$X4)))
F27alt <- (sum(F27$X4)/(sum(F27$X3) + sum(F27$X4)))
F50ref <- (sum(F50$X3)/(sum(F50$X3) + sum(F50$X4)))
F50alt <- (sum(F50$X4)/(sum(F50$X3) + sum(F50$X4)))
merge <- matrix(c(F27ref, F27alt, F50ref, F50alt), nrow=2, byrow=T)

fisher.test(merge)
# test rounds decimals to int so I'm just moving a decimal to test...
fisher.test((merge*10))
fisher.test((merge*100))


## Coverage
F27ref <- (sum(F27$X3)/nrow(F27))
F27alt <- (sum(F27$X4)/nrow(F27))
F50ref <- (sum(F50$X3)/nrow(F50))
F50alt <- (sum(F50$X4)/nrow(F50))
merge <- matrix(c(F27ref, F27alt, F50ref, F50alt), nrow=2, byrow=T)

fisher.test(merge)
fisher.test(merge*10)
fisher.test(merge*100)


t.test(F27$X3, F27$X4)
t.test(F50$X3, F50$X4)
t.test(F27$X3, F50$X3)
t.test((F27$X3+F27$X4), (F50$X3+F50$X4))

# allele counts relative to pop's total coverage
t.test(F27$X3/(F27ref), F50$X3/(F50ref))
t.test(F27$X4/(F27alt), F50$X4/(F50alt))

# combined pop counts relative to coverage
t.test((F27$X3+F27$X4)/(F27ref+F27alt), (F50$X3+F50$X4)/(F50ref+F50alt))
t.test((F27$X3+F27$X4)*10/(F27ref+F27alt), (F50$X3+F50$X4)*10/(F50ref+F50alt))




########## PLOTTING FRACTION OF SEGREGATING SITES OVER TIME ##########
# fraction of segregating sites over time for davis and bozeman
# include list of coverage for the runs
# plot previous and current segregating sites and coverages

## Previous run
generation_no <- c(0, 18, 28, 58, 18, 27, 50)
segregating_fraction_I <- c(1, 0.5867745,0.5232385,0.2770383, 0.5867745,0.4609283, 0.5630854)
location_I <- c("Davis", "Davis", "Davis","Davis", "Bozeman","Bozeman","Bozeman")
generation_name <- c("parents", "F18", "F28", "F58", "F18", "F27", "F50")


## BOZ II Sequencing
segregating_fraction_II <- c(1, 0.5867745,0.5232385,0.2770383, 0.5867745,0.12192, 0.1344696)
location_II <- c("Davis - 30x", "Davis - 30x", "Davis - 30x","Davis - 30x", "Bozeman - 3.5","Bozeman - 3.5","Bozeman - 3.5")
d <- tibble(generation_no, segregating_fraction_I, segregating_fraction_II, location_I, location_II, generation_name)

ggplot(d, aes(x=generation_no, y=segregating_fraction_I, group=location_I))+
  geom_line(aes(linetype=location_I))+
  geom_point(aes(shape=location_I))+
  theme_bw() 


ggplot(d, aes(x=generation_no, y=segregating_fraction_II, group=location_II))+
  geom_line(aes(linetype=location_II))+
  geom_point(aes(shape=location_II))+
  theme_bw()



