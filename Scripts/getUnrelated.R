
cat("Extracting unrelated IDs...\n")

fin="Files/20130606_g1k.ped"

sam=read.table(fin, stringsAsFac=F, head=T, sep="\t")

#pops to extract are
# LWK TSI CHB PEL
pops=c("LWK", "TSI", "CHB", "PEL")

cat("Written these files: ")
for (pop in pops) {
	pid=sam$Individual.ID[which(sam$Pop==pop & sam$Rela %in% c("mother","father","unrel","unrels"))]
	cat(pid, sep="\n", file=paste(pop, ".txt", sep="",collapse=""))
	#cat(pop, ":", length(pid), "\n")
	cat(paste(pop, ".txt", sep="",collapse=""), "")
}
cat("\n")



