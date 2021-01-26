setwd("/rhome/dkoenig/bigdata/BARLEY_CCII_POOLSEQ_WGS/DATA/OUTPUT/CALL_VARIANTS")
##############################################################################################
##############################################################################################
##############################################################################################
afs<-read.table("FULL_AFS.txt",stringsAsFactors = F)
#
pout<-rep(NA,nrow(afs))
for (i in 1:nrow(afs)){
  if (i/10000 == round(i/10000)){
    print (i)
  }
  pout[i]<-fisher.test(matrix(as.numeric(afs[i,5:12]),nrow=2, byrow = F))$p.value
}
afsout<-cbind(afs,pout)
write.table(afsout,"FULL_AFS_PVAL.txt",quote=F,row.names = F,col.names=F)
