
plotGenes<-function(chr,min.pos,max.pos){

  defaultPar<-par()

  xforrs = 0.03
  regsz = 1.2
  width=22

  rsqplus = 0.045
  rightylabxplus=0.05
  xmargin = 0.005
  cexaxis=1
  cexlab=1.5
  adj = 0

  filein="/home/matteo/Copenhagen/refGene.txt.gz"

  dat<-read.table(filein,as.is=T,head=T,comment.char="")
  xx2 = dat[dat[,"chrom"]==paste("chr",chr,sep="") & dat[,"cdsStart"]<max.pos*1e6 & dat[,"cdsEnd"] >min.pos*1e6,]

  start = xx2$txStart
  end   = xx2$txEnd
  nams  = xx2$name2
  cnts  = xx2$exonCount

  #par(mar=c(5.2, 6.2, -0.1, 6.3) +0.1)

  plot(c(0,0),c(0,0),type="n",xlim=c(min.pos-adj,max.pos+adj),ylim=c(-0.8,0.1),xlab="",
       xaxs="i",yaxt='n',ylab="",main="",cex.lab=1,cex.axis=cexaxis-0.4,tck=-0.05,
       mgp=c(3,0,0),las=1)
  mtext(1,text=paste("Position on chromosome ",chr," (Mb)",sep=""),line=2,cex=1)

  ord <- order(start)
  start    <- start[ord]
  end      <- end[ord]
  exoncnts <- cnts[ord]
  nams     <- nams[ord]
  keep <- !duplicated(nams)
  start    <- start[keep]
  end      <- end[keep]
  exoncnts <- cnts[keep]
  nams     <- nams[keep]
  ord <- ord[keep]
  he       <- rep(c(0,-0.18,-0.36,-0.54,-0.72),100)[1:length(nams)]-0.05

  if(length(start)>0){
    segments(start/1e6, he, end/1e6, he)
    keep = !duplicated(nams)
    sapply(1:sum(keep),function(x){text((end[keep][x]+start[keep][x])/2e6,he[keep][x]+0.08,bquote(italic(.(nams[keep][x]))),cex=cexlab-1.2)})
    estart = as.numeric(unlist(sapply(xx2$exonStarts[ord],function(y){strsplit(y,",")[[1]]})))/1e6
    eend = as.numeric(unlist(sapply(xx2$exonEnds[ord],function(y){strsplit(y,",")[[1]]})))/1e6
    rect(estart,rep(he,xx2$exonCount[ord])-0.01,eend, rep(he,xx2$exonCount[ord])+0.01,col="black")
  }

  par(mar=defaultPar$mar)
}


# adapted from a script by Javier Mendoza Revilla

args=commandArgs(T)
fin=args[1]
fout=args[2]
rm(args)

# source("Scripts/plotGenes.R")

pdf(file=fout)
par(mfrow=c(2,1))
par(mar=c(0, 4, 4, 2) + 0.1)

nsl <- read.table(fin, header=FALSE, stringsAsFact=F)

cat("Maximum nsL value:", max(nsl[,7], na.rm=T), "\n")
cat("Minimum nsL value:", min(nsl[,7], na.rm=T), "\n")

len=5e4
st=1e4

starts=seq(109e6,110e6-len,st)
ends=starts+len

vals=pos=rep(NA, length(starts))

for (i in 1:(length(starts))) {

	ind=which(nsl[,2]>=starts[i] & nsl[,2]<ends[i])
	vals[i]=length(which(abs(nsl[ind,7])>2))
	pos[i]=starts[i]+(len/2)

}


plot(pos, vals, type ="l", col="orange",ylab = "#|nSL|>2", xlab="Chromosome", main="nSL scan" , lwd=1, xaxt="n",cex.main = 0.9, cex.axis = 0.6, cex.lab = 0.68 ) #orange

#Add a second panel showing the genes in that region
par(mar=c(5, 4, 0.5, 2) + 0.1)
plotGenes(2,109,110)





