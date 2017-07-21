
# from the whole genetic map, extract sites in the VCF file, and produce file for input in selscan
# < chr # > <id > < genetic position > < physical position >
# whitespace delimited 

args=commandArgs(T)
fmap=args[1]
fvcf=args[2]
rm(args)


# fmap="Files/genetic_map_GRCh37_chr11.map"
# fvcf="Data/PEL.chr11.vcf"

vcf=read.table(fvcf, stringsAsFactors=F)
rs=vcf[,3]
pos=as.numeric(vcf[,2])
rm(vcf)

map=read.table(fmap, head=T)

ind=match(pos, map$Position)

maps=map$Map[ind]

app=approx(x=c(pos[1]-200,pos,pos[length(pos)]+200), y=c(min(maps,na.rm=T)-0.001, maps, max(maps,na.rm=T)+0.001), xout=c(pos[1]-200,pos,pos[length(pos)]+200))

newmap=cbind(rep(11, length(rs)), rs, app$y[2:(length(app$y)-1)], pos)

write.table(newmap, quote=F, sep=" ", col.names=F, row.names=F)



