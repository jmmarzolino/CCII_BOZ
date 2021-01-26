colnames(all_freqs_sub)<- c("chr", "position","0", "18", "27", "28", "50", "58")

test <- melt(all_freqs_sub,id=c("chr","position"))

ggplot(test,aes(x=position,y = value)) + geom_point()



# Lines and points; colour depends on cond2
ggplot(df2, aes(x=position, y=value)) + 
  geom_line(aes(colour=cond2, group=cond2)) + # colour, group both depend on cond2
  geom_point(aes(colour=cond2),               # colour depends on cond2
             size=3)                          # larger points, different shape
## Equivalent to above; but move "colour=cond2" into the global aes() mapping
# ggplot(df2, aes(x=cond1, y=yval, colour=cond2)) + 
#    geom_line(aes(group=cond2)) +
#    geom_point(size=3)