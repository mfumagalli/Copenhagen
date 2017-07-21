
## FROM MS FILE TO SUMMARY STATS

# not a general way to do it! it assumes you have 80 chroms and 4 pops with 20 chroms each!

args=commandArgs(T)
simfile=args[1]
rm(args)

source("/home/matteo/Copenhagen/Scripts/popgen.R")

nsam=80 # change this eventually
haplos=readMs2(simfile, nsam)

nrep=length(haplos$hap)

he=c("S","SS","TajimaD","FuLiDs","FuLiFs","H1", "H2", "H2/H1", "FST12", "FST13", "FST23", "PBS3")
cat(he, sep="\t")
cat("\n")
for (i in 1:nrep) {
	
		hap=haplos[[1]][[i]]

		haplist=list()
		haplist[[1]]=hap[21:40] # TSI
		haplist[[2]]=hap[41:60] # CHB
		haplist[[3]]=hap[61:80] # PEL
		rm(hap)
		# anc=paste(rep("0",nchar(haplist[[1]][1])),sep="",collapse="")

		## SUMMARY STATS
		taj=tajima(haplist[[3]])[5] #TD
		fu=fuli(haplist[[3]]) # 1 and 2 are S and SS, 3 and 4 are Ds and Fs
		hs=homohapl(haplist[[3]]) # H1, H2, H2/H1
		fsts=reynolds(haplist)
		pbs=dopbs(fsts[2],fsts[3],fsts[1])

		subres=c(fu[1:2],taj,fu[3:4],hs,fsts, pbs)
		cat(subres, sep="\t")
		cat("\n")
}





