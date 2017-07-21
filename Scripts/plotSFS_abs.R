
# from ANGSD website, adapted
# this is written knowning that the path is Results/

fins=commandArgs(T)

pops=unlist(strsplit(unlist(strsplit(fins, split="Results/")), split=".sfs"))
cat("Populations:", pops, "\n")
fout=paste("Results/", paste(pops, sep="",collapse="_"), ".pdf", sep="", collapse="")

pdf(file=fout)

msfs=c()
for (i in 1:length(fins)) {
	#read data
	sfs <- (scan(fins[i], quiet=T))
	cat(length(sfs)," ")
	#the variable categories of the sfs
	sfs <- sfs[-c(1,length(sfs))] 
	msfs <- rbind(msfs, sfs)
}

barplot(msfs, beside=T, legend=pops, xlab="Chromosomes", names=1:length(sfs), ylab="Counts", main="SFS")

dev.off()

cat("Output file:",fout,"\n")


