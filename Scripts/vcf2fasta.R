
# vcf to fasta

args=commandArgs(T)
fin=args[1]
pop=args[2]
fout=args[3] # this for site annotation
rm(args)

# append file to stdout

vcf=read.table(fin, sep="\t", stringsAsFact=F)
nind=ncol(vcf)-9
nsites=nrow(vcf)
vcf=vcf[,c(1:5,10:ncol(vcf))]

if (!is.null(fout)) {
	dim(vcf)
	write.table(vcf[,1:5], sep="\t", quote=F, row.names=F, col.names=F, file=fout)
}

## creata fasta
fas0=fas1=c()
for (i in 1:nind) fas0[i]=fas1[i]=""

for (i in 1:nsites) {
	        ref=vcf[i,4]; alt=vcf[i,5];
		for (j in 1:nind) {
			s=j+5
			all=strsplit(vcf[i,s],split="|")[[1]]
			if (as.numeric(all[1])==0) fas0[j]=paste(fas0[j], ref, sep="", collapse="") else fas0[j]=paste(fas0[j], alt, sep="", collapse="")
			if (as.numeric(all[3])==0) fas1[j]=paste(fas1[j], ref, sep="", collapse="") else fas1[j]=paste(fas1[j], alt, sep="", collapse="")
		}
}

for (i in 1:nind) {
	# pop
	cat(paste(">",pop,sep="",collapse=""),"\n")
        cat(fas0[i],"\n")
	cat(paste(">",pop,sep="",collapse=""),"\n")
	cat(fas1[i],"\n")
}


