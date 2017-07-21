

# function given by Javier Mendoza Revilla, possibly implemented by Emil Jørsboe inspired by LocusZoom

# gene annotation refer to hg19, downloaded http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/, 14th Feb 2016

# this function plots then genes for a certain interval on a chromosome
# it uses a hg19 gene list and some trick with par & plot
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

  filein="/home/matteo/Copenhagen/Files/refGene.txt.gz"

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

# paste above
# source("Scripts/plotGenes.R")

pdf(file=fout)
par(mfrow=c(2,1))
par(mar=c(0, 4, 4, 2) + 0.1)

pbs <- read.table(fin, header=TRUE)

pbs$PBS1[which(pbs$PBS1<0)]=0
pbs$PBS2[which(pbs$PBS2<0)]=0
pbs$PBS0[which(pbs$PBS0<0)]=0

cat("Maximum PBS value:", max(pbs$PBS2[which(pbs$midPos>109.4e6 & pbs$midPos<109.8e6)], na.rm=T), "\n")

plot(pbs$midPos/1e6, pbs$PBS2, type ="l", col="orange",ylab = "PBS", xlab=paste("Chromosome",pbs$chr[1]), main="PBS scan" , lwd=1, xaxt="n",cex.main = 0.9, cex.axis = 0.6, cex.lab = 0.68, ylim=c(0, max(pbs$PBS2, na.rm=T)), xlim=c(109.4, 109.8) ) #orange

lines(pbs$midPos/1e6, pbs$PBS0, col="blue", lwd=1)
#cat("Maximum PBS value (for LWK):", max(pbs$PBS0, na.rm=T), "\n")

legend("topright", legen=c("Target","Control"), col=c("orange", "blue"), lty=1)

    #title(ylab="Selection Statistics", line=2.2, cex.lab=0.69)
    #abline(h=1,lty=2,col="grey",lwd=1)


#Add a second panel showing the genes in that region
par(mar=c(5, 4, 0.5, 2) + 0.1)
plotGenes(pbs$chr[1],109.4,109.8)





