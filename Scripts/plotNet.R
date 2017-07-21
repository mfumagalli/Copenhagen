
library(ape)
library(ade4)
library(adegenet)
library(pegas)

args=commandArgs(T)
fin=args[1]
fin2=args[2]
fout=args[3]
rm(args)

x=read.dna(fin, format="fasta")

h=haplotype(x)
h=subset(h,1)
h <- sort(h, what = "labels")

countHap <- function(hap = h, dna = x){
	    with(
		         stack(setNames(attr(hap, "index"), rownames(hap))),
			         table(hap = ind, pop = attr(dna, "dimnames")[[1]][values])
			     )
}

ind.hap<-with(
	              stack(setNames(attr(h, "index"), rownames(h))),
		              table(hap=ind, pop=rownames(x)[values])
		              )

subset.haplotype <- function(x, freqmin = 1, freqmax = Inf, ...)
{
	    oc <- oldClass(x)
    idx <- attr(x, "index")
        f <- sapply(idx, length)
        s <- f <= freqmax & f >= freqmin
	    x <- x[s, ]
	    attr(x, "index") <- idx[s]
	        class(x) <- oc
	        x
}

net=haploNet(h)

pdf(file=fout)

plot(net, size=attr(net, "freq")*2, scale.ratio = 0.5, cex = 0.8, fast=TRUE, pie=ind.hap, labels=F, legend=F, show.mutation=0, th=0, lwd = (1 + round(.010*(net[,'step']))))
legend("topleft", colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20)

dev.off()

# identify mutations separating each haplotype

fas=readLines(fin)
fas=fas[seq(2,length(fas),2)]

anno=read.table(fin2, head=F)

ind=attr(h, "index")

for (i in 1:(length(ind)-1)) {
	for (j in (i+1):length(ind)) {
		s1=strsplit(fas[ind[[i]][1]],split="")[[1]]
		s2=strsplit(fas[ind[[j]][1]],split="")[[1]]
		idiff=anno[which(s1!=s2),2]
		cat(attr(h, "dimnames")[[1]][i], "\t", attr(h, "dimnames")[[1]][j], "\t", paste(idiff,sep="",collapse=","), "\n")
	}
}





