
args=commandArgs(T)

fin="Results/ALL.admix.K4.txt"

th=as.numeric(args[1])

pops=c(args[2], args[3], args[4], args[5])
inds=c(which(pops=="LWK"), which(pops=="TSI"), which(pops=="CHB"), which(pops=="NAM"))
rm(pops)

pops=c("AFR", "EUR", "EAS", "NAM")
cols=c("blue","red", "grey", "green")

pops=pops[inds]
cols=cols[inds]

pdf(file="Results/ALL.admix.pdf")

admix<-t(as.matrix(read.table(fin)[41:60,-1]))
barplot(admix,col=cols,space=0,border=NA,xlab="Individuals",ylab="admixture",legend=pops)

dev.off()

res=read.table(fin, stringsAsFactors=F, head=F)[41:60,]

# V1 is NAM, V2 is AFR, V3 is EUR

ii=which(res$V3>=th)
cat("Samples retained:", length(ii), "\n")

cat(res$V1[ii], sep="\n", file="Results/PEL.unadmixed.BAMs.txt")

cat("Output files:", "Results/ALL.admix.pdf", "Results/PEL.unadmixed.BAMs.txt", "\n")




