
# from ANGSD web site or from Thorfinn Korneliussen

getSFS <-function(x){
	a<-read.table(x)[,-c(1:2)]
 	ns=ncol(a)
	res=rep(0, (2*ns)+1  )
	a[a==-1] <- NA;
   	#tot=as.numeric(table(rowSums(a,na.rm=T)))
	tot=rowSums(a,na.rm=T)
	for (i in 1:length(res)) res[i]=length(which(tot==(i-1)))
     	#tot/sum(tot)
	res
}




fin=commandArgs(T)

cat(getSFS(fin))
cat("\n")


