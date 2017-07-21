
# plot histogram of simulated PBS and get p-value of observed data

args=commandArgs(T)
fin=args[1]
obs=as.numeric(args[2])
fout=args[3]
rm(args)

res=read.table(fin, head=T)
res$PBS3[res$PBS3<0]=0

pdf(file=fout)

hist(res$PBS3, breaks=30, xlab="PBS", main="Simulated values")
abline(v=obs, lty=2)

dev.off()

pv=(length(which(res$PBS3>=obs))+1)/nrow(res)
cat("P-value:", pv, "\n")





