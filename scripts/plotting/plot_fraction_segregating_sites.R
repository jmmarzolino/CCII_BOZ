library(tidyverse)
setwd("/bigdata/koeniglab/jmarz001/CCII_BOZ/results")

generation_no <- c(0, 18, 28, 58, 18, 27, 50)
segregating_fraction <- c(0.9997131, 0.5867745,0.5232385,0.2770383, 0.5867745,0.4609283, 0.5630854)
generation_name <- c("parents", "F18", "F28", "F58", "F18", "F27", "F50")
location <- c("Davis", "Davis", "Davis","Davis", "Bozeman","Bozeman","Bozeman")
d <- data.frame(generation_no, segregating_fraction, location, generation_name)


# plot davis and bozeman lines
ggplot(d, aes(x=generation_no, y=segregating_fraction, group=location, color=location)) + 
  geom_line(aes(), size=3) + 
  scale_color_manual(values=c("deepskyblue2", "darkorange1")) +
  geom_point(aes(shape=location), shape=c(16), size=5) + 
  theme_bw() +
  theme(text = element_text(size=20)) +
  ylab("Fraction of Sites Segregating") + 
  xlab("Generation") +
  ggtitle("Genetic Diversity in CC II")
ggsave("fraction_segregating_sites_DAV-BOZ.png", width = 9.5,height = 5, units = "in")



# plot just davis points and make it look nice
d %>%
  filter(location=="Davis") %>%
  ggplot(aes(x=generation_no, y=segregating_fraction, group=location, color=location)) + 
  geom_line(aes(), size=3) + 
  scale_color_manual(values="darkorange1") +
  geom_point(aes(shape=location), shape=c(16), size=5) + 
  theme_bw() +
  theme(text = element_text(size=20)) +
  ylab("Fraction of Sites Segregating") + 
  xlab("Generation") +
  ggtitle("Genetic Diversity in CC II")
ggsave("fraction_segregating_sites_DAV.png", width = 9.5,height = 5, units = "in")
  
  