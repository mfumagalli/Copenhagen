
args=commandArgs(T)
fin=args[1]
fout=args[2]
rm(args)

mds=read.table(fin, stringsAsFac=F, head=T)

pops<-c("LWK","TSI","PEL")
cols=c("red","green","blue")

print(cbind(pops,cols))

par(mfrow=c(1,2))

plot(mds[,1],mds[,2], col=rep(cols, each=20), pch=16, xlab="PC1", ylab="PC2")
plot(mds[,3],mds[,4], col=rep(cols, each=20), pch=16, xlab="PC3", ylab="PC4")



