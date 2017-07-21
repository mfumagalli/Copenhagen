
# random R functions
# use at your own risks!

dopbs<-function(Fgc,Fge,Fce) {
        Tgc= -log(1-Fgc)
        Tge= -log(1-Fge)
        Tce= -log(1-Fce)
        pbs= (Tgc + Tge - Tce)/2
        pbs
}



reynolds<-function(haplos, verbose=FALSE) {

 # calcola numerosità campione
 npop=length(haplos)
 new_nsam=c(); for (p in 1:npop) new_nsam[p]=length(haplos[[p]])

 # output
 reyn=c()

 # calcola frequenze alleliche
 sfs=countFreq(haplos, haplos[[1]][1], fixed.na=FALSE, plot=FALSE)

 t=0
 for (p1 in 1:(npop-1)) {
  for (p2 in (p1+1):npop) {

   t=t+1
   somma=0; sommaden=0
   n1=new_nsam[p1]; n2=new_nsam[p2]

   for (i in 1:length(sfs[[1]])) {

    pl1=sfs[[p1]][i]
    pl2=sfs[[p2]][i]
    if (verbose) cat("\n\npl",i,pl1,pl2)

    alfa1=1-((pl1^2)+((1-pl1)^2))
    alfa2=1-((pl2^2)+((1-pl2)^2))

    if (verbose) cat("\nalfa", alfa1, alfa2)

    Al = (0.5*(((pl1-pl2)^2)+(((1-pl1)-(1-pl2))^2))) - (((n1+n2)*(n1*alfa1+n2*alfa2)) / ((4*n1*n2)*(n1+n2-1)))

    AlBl= (0.5*(((pl1-pl2)^2)+(((1-pl1)-(1-pl2))^2))) + (((4*n1*n2 - n1 - n2)*(n1*alfa1 + n2*alfa2)) / ((4*n1*n2)*(n1+n2-1)))

    if (verbose) cat("\nA\tB\tA+B", Al, AlBl-Al, AlBl)

    if (!is.na(Al) & !is.na(AlBl)) {
     somma=somma+Al
     sommaden=sommaden+AlBl
    }

   }

# se 1 pop e' singleton e altra fissa e' 0
   if (somma==0 & sommaden==0) reyn[t]=NA else reyn[t]=somma/sommaden

  }
 }
#
#  cat(" NA:",length(which(is.na(reyn[g,]))))
#  cat(" ZERO:",length(which(reyn[g,]==0)))

 reyn

} # fine function


homohapl<-function(haplos) {

	# haplo homo
	H1=H12=H2H1=H2=NA
	hfreqs=sort(as.numeric(table(unlist(haplos)))/length(haplos), dec=T)
	if (length(hfreqs)>0) H1=sum(hfreqs^2)
	if (length(hfreqs)>1) {
        	H12=H1+(2*hfreqs[1]*hfreqs[2])
        	H2=H1-(hfreqs[1]^2)
        	H2H1=H2/H1
	}

	c(H1,H2,H2H1)

}


trim<-function (sequence, inner = FALSE) {
## Author: Giorgia Menozzi
    sequence <- sub("^ +", "", sequence)
    sequence <- sub(" +$", "", sequence)
    if (inner) {
        sequence <- gsub(" +", "", sequence)
    }
    sequence
}

readMs_sub<-function(filein, nsam, len, freg) {

 lfile<-readLines(con=filein)

 # start, end, pos indexes
 inds<-which(lfile=="//")+3
 inde<-which(lfile=="//")+3+nsam-1
 indpos<-inds-1

 # reg
 reg=read.table(freg, stringsAsFactors=F)
 valid=c()
 for (j in 1:nrow(reg)) valid=c(valid, seq(reg[j,1],reg[j,2]))
 valid=sort(unique(valid))

 res<-pos<-c()
 pos=round(as.numeric(strsplit(lfile[indpos],split=" ")[[1]][-1])*len)
 ind_valid=which(pos %in% valid)
 pos=pos[ind_valid]
 rm(reg)
 rm(valid)

 for (i in inds:inde) res<-c(res,paste(strsplit(lfile[i],split="")[[1]][ind_valid],sep="",collapse=""))

 readMs_sub<-list(hap=res, pos=pos)

}

readMs_pos<-function(filein, nsam, len) {

 lfile<-readLines(con=filein)

 # start, end, pos indexes
 inds<-which(lfile=="//")+3
 inde<-which(lfile=="//")+3+nsam-1
 indpos<-inds-1

 # reg
 pos=round(as.numeric(strsplit(lfile[indpos],split=" ")[[1]][-1])*len)

 res=c()
 for (i in inds:inde) res<-c(res,paste(strsplit(lfile[i],split="")[[1]],sep="",collapse=""))

 readMs_pos<-list(hap=res, pos=pos)

}


readMs2<-function(filein, nsam, len=c()) {

 lfile<-readLines(con=filein)
 inds<-which(lfile=="//")+3
 inde<-which(lfile=="//")+3+nsam-1
 # cat("Sample:",length(inds),"\n")
 indpos<-inds-1

 res<-pos<-list()
 for (i in 1:length(inds)) {
  res[[i]]<-trim(lfile[inds[i]:inde[i]])
  ss<-strsplit(lfile[indpos[i]],split=" ")[[1]]
  pos[[i]]<-as.numeric(ss[2:length(ss)])
  if (length(len)>0) pos[[i]]<-as.integer(pos[[i]]*len)
 }

 readMs2<-list(hap=res, pos=pos)

}


hap2ms<-function(haplos, outgroup, pos, fname, start, end) {

	if (is.list(haplos)) haplos=unlist(haplos)
	nchroms = length(haplos)
	len = nchar(haplos[1])
	outgroup = strsplit(outgroup, split="")[[1]]

	haps=c()
	for (i in 1:nchroms) {
		tmp=strsplit(haplos[i], split="")[[1]]
		newhap=rep(1,len)
		newhap[which(tmp==outgroup)]=0
		haps[i] =  paste(newhap, sep="",collapse="")
	}

	pos=(pos-start)/(end-start)
	pos=round(pos*(end-start))/(end-start)

	cat("-ms ", nchroms  ," 1 -s ", len, "\n0x2a6f029f03c2fcaf\n\n//\nsegsites: ", len, "\npositions: ", paste(pos, sep="", collapse=" "), "\n", paste(haps, sep="", collapse="\n"), sep="", collapse="", file=fname)

}



countFreq<-function(haplos, outgroup=c(), plot=F, fixed.na=FALSE) {

 if (!is.list(haplos)) haplos<-list(haplos)
 npop<-length(haplos)

 maf<-list()

 if (plot) { x11(); par(mfrow=c(npop,1)) }
 for (p in 1:npop) {

  nsam<-length(haplos[[p]])
  len<-nchar(haplos[[p]][1])
  maf[[p]]<-rep(NA, len)
  for (l in 1:len) {
   hs<-substring(haplos[[p]],l,l)
   hs<-hs[which(hs!="N")]
   if (length(outgroup)>0) {
    maf[[p]][l]<-length(which(hs!=substring(outgroup,l,l)))/length(which(hs!="N"))
   } else {
    maf[[p]][l]<-min(table(hs))/length(which(hs!="N"))
   }
  }
  if (fixed.na==TRUE) maf[[p]][which(maf[[p]]==1)]<-NA
  if (plot) if (length(outgroup)==0) hist(maf[[p]],breaks=10,xlim=c(0,0.5),sub=paste("nsam/2 is",nsam/2),main="MAF histogram") else hist(maf[[p]],breaks=10,xlim=c(0,1),sub=paste("nsam/2 is",nsam/2),main="DAF histogram")
 }

 countFreq<-maf

}

countDaf<-function(haplos, outgroup, fixed.na=FALSE) {

 if (!is.list(haplos)) haplos<-list(haplos)
 npop<-length(haplos)

 maf<-list()

 for (p in 1:npop) {

  nsam<-length(haplos[[p]])
  len<-nchar(haplos[[p]][1])
  maf[[p]]<-rep(NA, len)
  for (l in 1:len) {
   if (substring(outgroup,l,l)!="N") {
    hs<-substring(haplos[[p]],l,l)
    maf[[p]][l]<-length(which(hs!=substring(outgroup,l,l)))/length(which(hs!="N"))
    }
  }
  if (fixed.na==TRUE) maf[[p]][which(maf[[p]]==1)]<-NA
 }

 countFreq<-maf

}




faywuH<-function(daf, nsam) {

  # daf is array of derived alle freq, for one single pop, per each site, output of countFreq

  thetaH=thetaP=FWU=0;

  for (i in 1:(nsam-1)) {

    Si=length(which(daf==i))
    thetaH<-thetaH + ( (2*Si*(i^2)) / (nsam*(nsam-1)) )
    thetaP<-thetaP + ( (2*Si*i*(nsam-i)) / (nsam*(nsam-1)) )
  }

  FWH<-thetaP-thetaH
 
  FWH
}

colorSweep<-function(haplos, outgroup=NA, mut, show.labels=TRUE, show.grid=TRUE, show.ticks=TRUE, main="") {

 if (!is.list(haplos)) haplos<-list(haplos)

 npop<-length(haplos)
 nsam<-rep(NA,npop); for (p in 1:npop) nsam[p]<-length(haplos[[p]])
 csam<-c(0, cumsum(nsam))
 ipop<-list(); for (p in 1:(length(csam)-1)) ipop[[p]]<-seq(csam[p]+1,csam[p+1])

 len<-nchar(haplos[[1]][1]) # nr siti

 lanci<-lder<-c()

 # trasforma aplotipi in formato ms (0 e 1 in relazione all'ancestore)
 if (!is.na(outgroup)) {
  f2m<-fas2ms(haplos=unlist(haplos), outgroup=outgroup, snp_pos=seq(1,nchar(unlist(haplos)[1])))
  apli<-f2m[[1]][[1]]
  muti<-f2m[[2]][[1]]
 } else {
  apli<-haplos
  muti<-seq(1,len)
 }

 mut<-which(muti==mut)
 cat("mut:",mut)

 anci<-which(substring(apli,mut,mut)=="0")
 der<-which(substring(apli,mut,mut)=="1")
 cat("\tanc:",length(anci)," ","der:",length(der))

 # per ogni pop
 nanc<-nder<-rep(NA, npop)
 for (p in 1:npop) {
  nanc[p]<-length(which(!is.na(match(anci, ipop[[p]]))))
  nder[p]<-length(which(!is.na(match(der, ipop[[p]]))))
 }
 cat("\nnanc:",nanc,"\nnder:",nder) 

 lanci<-strsplit(apli[anci],split="")
 lder<-strsplit(apli[der],split="")

 lapl<-c(lanci,lder)

 # crea input per image
# xx<-seq(0,len)
 xx<-seq(0,nchar(apli[1]))
 yy<-seq(0,length(unlist(haplos)))
# zz<-matrix(0,ncol=length(lapl),nrow=len)
 zz<-matrix(0,ncol=length(lapl),nrow=nchar(apli[1]))
 for (i in 1:length(lapl)) {
  ind1<-which(lapl[[i]]=="1")
  zz[ind1,i]<-1
 }
 im<-list(x=xx,y=yy,z=zz)

 # plot
 #x11()
 if (show.ticks) par(lab=c(len,csam[length(csam)],4)) else par(xaxt="n", yaxt="n")
 if (show.labels) {
#  par(lab=c(len,csam[length(csam)],4))
  image(im, xlab="Segregating site", ylab="Sample", main=main, sub=paste("mut:",mut,"anc:",length(anci),"der:",length(der),"; haplos anc:", paste(nanc,sep=" ", collapse=" "),"der",paste(nder,sep=" ",collapse=" ")))
 } else {
#  par(lab=c(len,csam[length(csam)],0))
  image(im, main=main)
 }
 if (show.grid) {
  for (i in 1:csam[length(csam)]) lines(c(0,len),c(i,i), lty=1, lwd=0.3)
  for (i in 1:len) lines(c(i,i),c(0,csam[length(csam)]), lty=1, lwd=0.3)
  lines(c(0,len),c(length(anci),length(anci)), lty=1, lwd=4)
 }

} # fine function


dind<-function(haplos, outgroup, snp_pos=c(), addInteger=20) {

 # ref. barreiro et al. plos genetics 2009, toll-like receptor

 if (length(snp_pos)==0) snp_pos<-seq(1,nchar(haplos[1]))

 # to get ancestral/derived state anc-der
 newhap<-fas2ms(haplos=haplos, outgroup=outgroup, snp_pos=snp_pos)
  newsnp<-newhap[[2]][[1]]
  newhap<-newhap[[1]][[1]]

# daf<-countFreq(haplos, outgroup, plot=F, fixed=T)[[1]]
# daf<-daf[which(!is.na(daf))]

 Dind<-daf<-Dind_anc<-rep(NA, nchar(newhap[1]))

 # x ogni snp calcola rapporto tra k der e k anc
 for (s in 1:nchar(newhap[1])) {

  iA<-which(substring(newhap,s,s)=="0")
  iD<-which(substring(newhap,s,s)=="1")

  daf[s]<-length(iD)/length(haplos)

 # K (nucleotide differences)
  Ahap<-newhap[iA]
  n<-length(iA)
  if (n==1) AKl<-0 else { 
  Akij<-0
   for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      Akij<-Akij + length(which(strsplit(Ahap[i],split="")[[1]]!=strsplit(Ahap[j],split="")[[1]]))
    }
   }
   AKl<-Akij/choose(n,2)
  }

 # K (nucleotide differences)
  Dhap<-newhap[iD]
  n<-length(iD)
  if (n==1) DKl<-0 else {
   Dkij<-0
   for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      Dkij<-Dkij + length(which(strsplit(Dhap[i],split="")[[1]]!=strsplit(Dhap[j],split="")[[1]]))
    }
   }
   DKl<-Dkij/choose(n,2)
  }

 # dind
 if (DKl!=0) Dind[s]<-AKl/DKl
 if (AKl!=0) Dind_anc[s]<-DKl/AKl

 }

 Dind[which(is.na(Dind))]<-max(Dind, na.rm=T)+addInteger
 Dind_anc[which(is.na(Dind_anc))]<-max(Dind_anc, na.rm=T)+addInteger
 
# if (is.na(addInteger)) {
#   Dind[which(is.na(Dind))]<-0
#   Dind_anc[which(is.na(Dind_anc))]<-0
# }

 res<-list(snp_pos=newsnp, daf=daf, dind=Dind, dind_anc=Dind_anc)

 res

}

fas2ms<-function(haplos, outgroup, fname="", snp_pos=NA, lenseq=NA, verbose=F) {

 if (!is.list(haplos)) haplos<-list(haplos)
 npop<-length(haplos)

 res<-list()
 snps<-list()

 for (p in 1:npop) {

  snps_tmp<-c()

  nsam<-length(haplos[[p]])
  len<-nchar(haplos[[p]][1])
  newhaplos<-rep("",nsam)

  notvalid<-0

  for (l in 1:len) {

   alleli<-substring(haplos[[p]],l,l)
   anc<-substring(outgroup,l,l)

   if ( (length(unique(alleli))==2) & !is.na(match(anc,unique(alleli))) & anc!="N" & anc!="-" )  { # se anc è in alleli

    snps_tmp<-c(snps_tmp, snp_pos[l])

    for (n in 1:nsam) {

     if(alleli[n]==anc) {
      newhaplos[n]<-paste(newhaplos[n],0,sep="")
     } else {
      newhaplos[n]<-paste(newhaplos[n],1,sep="")
     }

    }

   } else notvalid<-notvalid+1 #cat(" ",l) # se nn sono biallelici oppure non c'è anc

  }

  res[[p]]<-newhaplos
  snps[[p]]<-snps_tmp

  if (verbose) cat("Not valid:", notvalid, "\n")

 }

 if (nchar(fname)>0 & !is.na(snp_pos[1]) & npop==1) {
  snippi<-(snps[[1]]/lenseq)
  cat("",file=fname)
  cat("//\nsegsites:", length(snps[[1]]), "\npositions:",snippi,"\n",file=fname)
  cat(unlist(res[[1]]),sep="\n",file=fname,append=TRUE)
 }

 fas2ms<-list(res,snps)

}

linkage<-function(haplos, na.rm=FALSE, plot=F) {

 if (!is.list(haplos)) haplos<-list(haplos)
 npop<-length(haplos)

 res<-list()

 for (p in 1:npop) {

  if (na.rm) {
   polyss<-polyndex(haplos[[p]])$snp_pos
   for (i in 1:length(haplos[[p]])) {
    haplos[[p]][i]<-paste(strsplit(haplos[[p]][i],split="")[[1]][polyss],sep="",collapse="")
   }
  } else polyss<-NA

  len<-nchar(haplos[[p]][1])
  nsam<-length(haplos[[p]])
  cat("Sample:",nsam,"Sites:",len,"\n")
 # cat("\n",rep("",(floor(len/10)-1)),"|\n")
  Dprime<-rsquare<-pv<-Dprimeplot<-matrix(NA, nrow=len, ncol=len)

  for (i in 1:(len-1)) {

   for (j in (i+1):len) {

    # alleli
    subi<-substring(haplos[[p]],i,i)
    subj<-substring(haplos[[p]],j,j)
    # forme alleliche
    ali<-unique(subi)
    alj<-unique(subj)
    if (length(ali)==1) ali<-c(ali, "N") # se monomorfico, soluzione temporanea
    if (length(alj)==1) alj<-c(alj, "N")

    # freq alleliche
    pp<-rep(NA,4)
    pp[1]<-length(which(subi==ali[1]))/nsam
    pp[2]<-length(which(subi==ali[2]))/nsam
    pp[3]<-length(which(subj==alj[1]))/nsam
    pp[4]<-length(which(subj==alj[2]))/nsam

    # aplotipi
    haps<-paste(subi,subj,sep="")

    # poss aplotipi
    posshap<-c()
    for (ai in 1:2) {
     for (aj in 1:2) {
      posshap<-c(posshap, paste(ali[ai],alj[aj],sep="",collapse=""))
     }
    }

    # freq aplotipiche
    P<-rep(NA,4) # prob aplotipi
    for (c in 1:4) {
     P[c]<-length(which(haps==posshap[c]))/nsam
    }

    # D, Dprime
    D<-P[1]*P[4]-P[2]*P[3]
    Dmax<-min(c(pp[1]*pp[4]),c(pp[2]*pp[3]))
    Dmin<-max(c(-pp[1]*pp[3]),c(-pp[2]*pp[4]))
    if (D<0) {
     Dprime[i,j]<-abs(D/Dmin)
    } else {
     Dprime[i,j]<-abs(D/Dmax)
    }

    # rsquare
    rsquare[i,j]<-(D^2)/prod(pp)

    # test
    chiquadro<-rsquare[i,j]*nsam
    pv[i,j]<-pchisq(q=chiquadro, df=1, lower.tail=F)

   } # fine for in j
  } # fine for in i

  Dprimeplot[which(Dprime==1 & pv<=0.05)]<-0.1
  Dprimeplot[which(Dprime==1 & pv>0.05)]<-0.3
  Dprimeplot[which(Dprime<1 & pv<=0.05)]<-0.6
  Dprimeplot[which(Dprime<1 & pv>0.05)]<-0.9

 if (plot) {
  x11()
  par(mfrow=c(2,2))
  image(x=seq(1,len), y=seq(1,len),z=-(matrix(findInterval(Dprime,seq(0,1,0.2)),nrow=len,ncol=len)), col=heat.colors(5), main="Dprime", xlab="",ylab="", sub=paste("Pop",p))
  image(x=seq(1,len), y=seq(1,len),z=-(matrix(findInterval(rsquare,seq(0,1,0.2)),nrow=len,ncol=len)), col=heat.colors(5), main="R^2", xlab="",ylab="", sub=paste("Pop",p))
  image(x=seq(1,len), y=seq(1,len),z=(matrix(findInterval(pv,c(0,0.01,0.05,1)),nrow=len,ncol=len)), col=heat.colors(3), main="pvalue", xlab="",ylab="", sub=paste("Pop",p))
  image(x=seq(1,len), y=seq(1,len),z=(matrix(findInterval(Dprimeplot,seq(0,1,0.25)),nrow=len,ncol=len)), col=c("red","lightblue","lightpink","white"), main="Combined", xlab="",ylab="", sub=paste("Pop",p))

 }

  res[[p]]<-list(Dprime=Dprime, rsquare=rsquare, pv=pv, snp=polyss)

 } # fine for in p

 linkage<-res

} # fine function



EHHS<-function(haplos_tot, th=0.10, snp_pos=c(), index=c(), rsb=TRUE) {

 # input parameters:
 # haplos_tot: list of haplotypes
 # th: threshold for the decay of signal when computing integrating scores
 # snp_pos: position of polymorphic sites (you need it for integrating EHHS)
 # index: which SNPs are core haplotypes?
 # rsb: boolean, compute rsb?

 # haplotypes are as follow:
 # a list, each field of the list is an array of phased chromosomes, where 1st and 2nd are the two haplotypes for 1st individual, 3rd and 4th are the two haplotypes for 2nd individual and so on.
 # e.g.:
#  haplos=list()
#  haplos[[1]]=c("AGT","AGG","ATT","AGT")
#  haplos[[2]]=c("GGT","ATG","ATG","ATG")
 # 2 populations, 2 indivuals each, with 3 SNPs

 # for each site and each population compute iES through EHHS
 if (!is.list(haplos_tot)) haplos=list(haplos_tot)
 npop=length(haplos_tot)
 nsite=nchar(haplos_tot[[1]][1])

 # if not specified, test all SNPs
 if (length(index)==0) index=1:nsite

 # I need snp_pos to integrate EHHS
 if (length(snp_pos)==0) snp_pos=seq(1,nsite,1)

 # initialize
 ehhs=ehha=list()
 ies=matrix(0, nrow=npop, ncol=nsite)

 # cycle for each pop
 for (p in 1:npop) {

  ehhs[[p]]=matrix(0, nrow=nsite, ncol=nsite)
  ehha[[p]]=list()

  # cycle for each tested SNP
  for (i in index) {

   #cat("\n",p,", ",i,"/",nsite,":")

   haplos=substring(haplos_tot[[p]],i,i)

   alleli=unique(haplos)
   
   pa=c() # allele frequencies
 
   value_site=matrix(0, nrow=length(alleli), ncol=nsite)

   # going forward...
   for (s in 1:length(alleli)) {

    ehha[[p]][[s]]=matrix(0, nrow=nsite, ncol=nsite)

    icore=which(haplos==alleli[s]) # who has got 'core', thus first allele [chi ha core, ovvero primo allele]
    pa=c(pa,length(icore))    
    # denominator formula EHH
    denom<-choose(length(icore),2)
#     lcore<-1
#     lhaplo<-nchar(haplos)[1]
 #   haplos<-haplos[icore] # seleziono solo quelli che hanno il core
    # x ogni finestra calcolare il valore di EHH, step=1

    if (denom==0) {

     value_site[s,]=rep(0,ncol(value_site))

    } else {

    estremo<-seq(i,nsite,1)
    value<-rep(0,length(estremo))
    for (w in 1:length(estremo)) {
     hapi<-substring(haplos_tot[[p]][icore],i,estremo[w])
     ss<-length(unique(hapi)) # see sweep documentation for notation
     et<-as.numeric(table(hapi))
     numer<-0
     for (ii in 1:ss) {
      numer<-numer+choose(et[ii],2)
     }
     value[w]<-numer/denom
    } # end 'for' in w
    value_site[s,i:nsite]=value

    # backward
    estremo<-seq(i,1,-1)
    value<-rep(0,length(estremo))
    for (w in 1:length(estremo)) {
     hapi<-substring(haplos_tot[[p]][icore],estremo[w],i)
     ss<-length(unique(hapi)) # see sweep documentation for notation
     et<-as.numeric(table(hapi))
     numer<-0
     for (ii in 1:ss) {
      numer<-numer+choose(et[ii],2)
     }
     value[w]<-numer/denom
    } # end 'for' in w
    value_site[s,1:i]=rev(value)

   }

    ehha[[p]][[s]][i,]=value_site[s,]
    names(ehha[[p]])[s]=alleli[s]

   } # end 'for' s in alleli

   pa=pa/(sum(pa))
   # these are EHHA values, compute EHHS
   value_ehhs=apply(X=value_site, FUN=weighted.mean, MAR=2, w=pa^2)

   ehhs[[p]][i,]=value_ehhs

   # compute iES, reference tang 2007
   # retrieve where/when EHHS falls below our threshold (0.10)
   a=max(setdiff(which(value_ehhs[1:i]<th),i)); if (abs(a)==Inf) a=i
   b=min(setdiff(which(value_ehhs[i:nsite]<th),i))+i-1; if (abs(b)==Inf) b=i
   #cat(" a:",a,"b:",b)

   ties=0 # ies
   for (j in (a+1):b) {
    ties=ties+(value_ehhs[j-1]+value_ehhs[j])*((snp_pos[j]-snp_pos[j-1]))/2
   }

   ies[p,i]=(ties)
 
   #cat(" ",ties)

  } # end 'for' i in nsite

 } # end 'for' p in npop

 lnrsb=c()
 # compute rsb
 if (rsb==TRUE & npop>1) {

  for (p1 in 1:(npop-1)) {

   for (p2 in (p1+1):npop) {

    lnrsb=rbind(lnrsb, log(ies[p1,]/ies[p2,]))    
   
   } # end 'for' in p2

  } # end 'for' in p1

 } # end if compute differential

 res=list(ehha=ehha, ehhs=ehhs, ies=ies, lnrsb=lnrsb)

 res

} # end function




linkage<-function(haplos, na.rm=FALSE, plot=F) {

 if (!is.list(haplos)) haplos<-list(haplos)
 npop<-length(haplos)

 res<-list()

 for (p in 1:npop) {

  if (na.rm) {
   polyss<-polyndex(haplos[[p]])$snp_pos
   for (i in 1:length(haplos[[p]])) {
    haplos[[p]][i]<-paste(strsplit(haplos[[p]][i],split="")[[1]][polyss],sep="",collapse="")
   }
  } else polyss<-NA

  len<-nchar(haplos[[p]][1])
  nsam<-length(haplos[[p]])
  cat("Sample:",nsam,"Sites:",len,"\n")
 # cat("\n",rep("",(floor(len/10)-1)),"|\n")
  Dprime<-rsquare<-pv<-Dprimeplot<-matrix(NA, nrow=len, ncol=len)

  for (i in 1:(len-1)) {

   for (j in (i+1):len) {

    # alleli
    subi<-substring(haplos[[p]],i,i)
    subj<-substring(haplos[[p]],j,j)
    # forme alleliche
    ali<-unique(subi)
    alj<-unique(subj)
    if (length(ali)==1) ali<-c(ali, "N") # se monomorfico, soluzione temporanea
    if (length(alj)==1) alj<-c(alj, "N")

    # freq alleliche
    pp<-rep(NA,4)
    pp[1]<-length(which(subi==ali[1]))/nsam
    pp[2]<-length(which(subi==ali[2]))/nsam
    pp[3]<-length(which(subj==alj[1]))/nsam
    pp[4]<-length(which(subj==alj[2]))/nsam

    # aplotipi
    haps<-paste(subi,subj,sep="")

    # poss aplotipi
    posshap<-c()
    for (ai in 1:2) {
     for (aj in 1:2) {
      posshap<-c(posshap, paste(ali[ai],alj[aj],sep="",collapse=""))
     }
    }

    # freq aplotipiche
    P<-rep(NA,4) # prob aplotipi
    for (c in 1:4) {
     P[c]<-length(which(haps==posshap[c]))/nsam
    }

    # D, Dprime
    D<-P[1]*P[4]-P[2]*P[3]
    Dmax<-min(c(pp[1]*pp[4]),c(pp[2]*pp[3]))
    Dmin<-max(c(-pp[1]*pp[3]),c(-pp[2]*pp[4]))
    if (D<0) {
     Dprime[i,j]<-abs(D/Dmin)
    } else {
     Dprime[i,j]<-abs(D/Dmax)
    }

    # rsquare
    rsquare[i,j]<-(D^2)/prod(pp)

    # test
    chiquadro<-rsquare[i,j]*nsam
    pv[i,j]<-pchisq(q=chiquadro, df=1, lower.tail=F)

   } # fine for in j
  } # fine for in i

  Dprimeplot[which(Dprime==1 & pv<=0.05)]<-0.1
  Dprimeplot[which(Dprime==1 & pv>0.05)]<-0.3
  Dprimeplot[which(Dprime<1 & pv<=0.05)]<-0.6
  Dprimeplot[which(Dprime<1 & pv>0.05)]<-0.9

 if (plot) {
  x11()
  par(mfrow=c(2,2))
  image(x=seq(1,len), y=seq(1,len),z=-(matrix(findInterval(Dprime,seq(0,1,0.2)),nrow=len,ncol=len)), col=heat.colors(5), main="Dprime", xlab="",ylab="", sub=paste("Pop",p))
  image(x=seq(1,len), y=seq(1,len),z=-(matrix(findInterval(rsquare,seq(0,1,0.2)),nrow=len,ncol=len)), col=heat.colors(5), main="R^2", xlab="",ylab="", sub=paste("Pop",p))
  image(x=seq(1,len), y=seq(1,len),z=(matrix(findInterval(pv,c(0,0.01,0.05,1)),nrow=len,ncol=len)), col=heat.colors(3), main="pvalue", xlab="",ylab="", sub=paste("Pop",p))
  image(x=seq(1,len), y=seq(1,len),z=(matrix(findInterval(Dprimeplot,seq(0,1,0.25)),nrow=len,ncol=len)), col=c("red","lightblue","lightpink","white"), main="Combined", xlab="",ylab="", sub=paste("Pop",p))

 }

  res[[p]]<-list(Dprime=Dprime, rsquare=rsquare, pv=pv, snp=polyss)

 } # fine for in p

 linkage<-res

} # fine function

polyslow<-function(haplos) {

 # calcola le statistiche di neutralità
 # HELP
 # haplos: lista con aplotipi divisi per popolazione (get_haplo$haplist)
 # S: nr snp (get_haplo$ssite)
 # pos: pos snp (get_haplo$ssite_pos)

 # convertire in lista
 if (!is.list(haplos)) haplos<-list(haplos)

 ## per ogni popolazione
 npop<-length(haplos)

 # risultato è una matrice, stampa i valori per ogni opolazione/riga
 res<-matrix(NA, nrow=npop, ncol=11)
 colnames(res)<-c("nsam","len","S","Ss","Ws","Wl","Ks","Kl","tajimaD","fuliDs","fuliFs")

 for (p in 1:npop) {

 # cat("\nPop:",p)
  haplosi<-haplos[[p]] # aplotipi dell p-esima popolazione

  # polyndex
  pindex<-polyndex(haplos=haplosi)
  pos<-pindex[[2]]
  res[p,1:8]<-pindex[[1]]

  # tajima
  res[p,9]<-tajima(haplos=haplosi, S=res[p,3], Kl=res[p,8])[[4]]

  # fuli
  res[p,10:11]<-fuli(haplos=haplosi, S=res[p,3], Ss=res[p,4], Kl=res[p,8], pos=pos)[3:4]

 } # fine for in p

 res
 #cat("\n")
 #print(res)

}



polyndex <- function (haplos, S=NA, pos=NA, popname=NA) {

 # calcola gli indici e i coefficienti per le statistiche di neutralità a partire da una lista contenente gli aplotipi

 # riceve gli aplotipi
 # se si conoscono già (ad esempio dopo funzione get_haplo) immettere anche nr snp (S) e loro posizione relativa alla sequenza (pos)

  haplos<-toupper(haplos)

  # numerosità campioni = numero aplotipi
  nsam<-length(haplos)
  # lunghezza campione
  len<-nchar(haplos)[1]

  # 1) S: nr siti polimorfici
  if (is.na(S[1]) | is.na(pos[1])) {
   # calcola numero di segr.sites S
   numall<-rep(NA, len)
   for (l in 1:len) { # x ogni base
    numall[l]<-length(unique(substring(haplos,l,l)))
   } # fine for in i 
   pos<-which(numall>1) # posizioni con snp
   S<-length(pos)
  } # fine calcola S e pos

  # 2) Ss: nr singleton
  Ss<-0
  for (l in pos) {
   df<-as.data.frame(table(substring(haplos,l,l)))
   Ss<-Ss+length(which(df$Freq==1))
  }

  # 3) an bn (a1 a2 x tajima), Ws Wl (theta watterson per site or on total length) 
  an<-bn<-0
  for (n in 1:(nsam-1)) {
   an<-an+(1/n)
   bn<-bn+(1/(n^2))
  }
  Wl<-S/an
  Ws<-Wl/len

  # 4) K (nucleotide differences)
  kij<-0
  for (i in 1:(n-1)) {
   for (j in (i+1):n) {
     kij<-kij + length(which(strsplit(haplos[i],split="")[[1]]!=strsplit(haplos[j],split="")[[1]]))
   }
  }
  Kl<-kij/choose(n,2)
  Ks<-Kl/len

  # scrivi output
#  polyndex<-c(p, nsam, len, S, Ss, Ws, Wl, Ks, Kl)
 ris<-c(nsam, len, S, Ss, Ws, Wl, Ks, Kl)

 polyndex<-list(ris,pos)

} # fine function


tajima<-function(haplos, S=NA, Kl=NA) {

 # riceve gli aplotipi e calcola tajima's D (Tajima 1989)

 # numerosità campioni = numero aplotipi
 n<-length(haplos)
 # lunghezza campione
 len<-nchar(haplos)[1]

 # calcola numero di segr.sites S
 haplos<-toupper(haplos)

 # S
 if (is.na(S)) {
  numall<-rep(NA, len)
  for (l in 1:len) { # x ogni base
   numall[l]<-length(unique(substring(haplos,l,l)))
  } # fine for in i 
  pos<-which(numall>1) # posizioni con snp
  S<-length(pos)
 }

 # a1 a2
 a1<-a2<-0
  for (i in 1:(n-1)) {
  a1<-a1+(1/i)
  a2<-a2+(1/(i^2))
 }

 # b1 b2
 b1<-(n+1)/(3*(n-1))
 b2<-(2*((n^2)+n+3))/(9*n*(n-1))

 # c1 c2
 c1<-b1-(1/a1)
 c2<-b2-(n+2)/(a1*n)+a2/(a1^2)

 # e1 e2
 e1<-c1/a1
 e2<-c2/((a1^2)+a2)

 # k, nucleotide differences
 if (is.na(Kl[1])) {
  kij<-0
  for (i in 1:(n-1)) {
   for (j in (i+1):n) {
     kij<-kij + length(which(strsplit(haplos[i],split="")[[1]]!=strsplit(haplos[j],split="")[[1]]))
   }
  }
  Kl<-kij/choose(n,2)
 }

 # Tajima's D
 D <- (Kl-(S/a1)) / sqrt( e1*S + e2*S*(S-1) )

 # risultato
 #tajima<-list(pos=pos, S=S, W=S/a1, Kl=Kl, D=D)
 tajima<-c(0,S,S/a1,Kl,D)

} # fine function


fuli<-function(haplos, S=NA, Ss=NA, pos=NA, Kl=NA) {

 # calcola Fu Li's D*, F* (1993) con correzione per F* di Simonsen (1995).
 # riceve in ingresso gli aplotipi

 # numerosità campioni = numero aplotipi
 n<-length(haplos)

 # lunghezza campione
 len<-nchar(haplos)[1]

 # calcola numero di segr.sites S
 haplos<-toupper(haplos)

 if (is.na(S[1]) | is.na(pos[1])) {
  numall<-rep(NA, len)
  for (l in 1:len) { # x ogni base
   numall[l]<-length(unique(substring(haplos[1:n],l,l)))
  } # fine for in i 
  pos<-which(numall>1) # posizioni con snp
  S<-length(pos)
 }

 # Ss nr singleton
 if (is.na(Ss[1])) {
  Ss<-0
  for (l in pos) {
   df<-as.data.frame(table(substring(haplos[1:n],l,l)))
   Ss<-Ss+length(which(df$Freq==1))
  }
 }

 # an bn an1
 an<-bn<-an1<-0
 for (i in 1:(n-1)) {
  an<-an+(1/i)
  bn<-bn+(1/(i^2))
 }
 an1<-an+(1/n)

 # cn dn
 cn <- 2 * ( ((n*an)-2*(n-1)) / ((n-1)*(n-2)) )
 dn <- cn + ((n-2)/((n-1)^2)) + (2/(n-1)) * (3/2 - (2*an1-3)/(n-2) - 1/n) 

 # vds (v D star)
 vds <- ( ((n/(n-1))^2)*bn + (an^2)*dn - 2*(n*an*(an+1))/((n-1)^2) ) / (an^2+bn)
 # uds
 uds <- ( (n/(n-1)) * (an - n/(n-1)) ) - vds

 # Fu Li's Ds (Dstar)
 Ds <- ( (n/(n-1))*S - (an*Ss) ) / sqrt(uds*S + vds*(S^2))

 if (is.na(Kl)) {
  # k, nucleotide differences
  kij<-0
  for (i in 1:(n-1)) {
   for (j in (i+1):n) {
     kij<-kij + length(which(strsplit(haplos[i],split="")[[1]]!=strsplit(haplos[j],split="")[[1]]))
   }
  }
  Kl<-kij/choose(n,2)
 }

 # vfs simonsen 1995
 vfs <- ( (( 2*(n^3) + 110*(n^2) -255*n + 153 ) / ( 9*(n^2)*(n-1) )) + ((2*(n-1)*an)/(n^2)) - ((8*bn)/n) ) / ((an^2) + bn)

 # ufs
 ufs <- ( ( n/(n+1) + (n+1)/(3*(n-1)) - 4/(n*(n-1)) + ((2*(n+1))/((n-1)^2))*(an1-((2*n)/(n+1))) ) / an ) - vfs

 # Fu Li's Fstar
 Fs <- ( Kl - (((n-1)/n)*Ss) ) / sqrt(ufs*S + vfs*(S^2))

 # Fu Li
 fuli<-c(S, Ss, Ds, Fs)

} # fine function


fusfs<-function(haplos, stirlingFirstKind=NA) {

 if (!is.list(haplos)) haplos<-list(haplos)
 nmax<-0
 for (p in 1:length(haplos)) {
  nmax<-max(c(nmax, length(haplos[[p]])))
 }

 if (is.na(stirlingFirstKind[1])) {
  source("/adat/Software/R/rtx/stirling.rtx")
  stirlingFirstKind<-StirlingMat(nmax+1)
 }

 res<-rep(NA, length(haplos))

 cat(".")
 for (p in 1:length(haplos)) {

  # compute Fu's Fs statistics Fu 1997
  if ((p %% 1000)==0) cat("\n",p)
  pstat<-polyndex(haplos=haplos[[p]])$ris

  nsam <- as.numeric(pstat[1]) # nr sequences
  k0 <- length(unique(haplos[[p]])) # "allele" è nr aplotipi (unici)
  thetapi <- as.numeric(pstat[8]) # mean nr nucl pairwise diff

  denom<-1
  for (n in 0:(nsam-1)) {
   denom<-denom*(thetapi+n)
  }

  Sprime<-0
  for (k in k0:(nsam)) {
   Sprime<-Sprime + (abs(stirlingFirstKind[nsam+1,k+1])*(thetapi^k)/denom)
  }

  res[p]<-log(Sprime/(1-Sprime))

 } # fine for in p

 res

}

kellyZns<-function(haplos) {

 if (!is.list(haplos)) haplos<-list(haplos)

 npop<-length(haplos)

 res<-rep(NA, npop)

 for (p in 1:npop) {

  ll<-linkage(haplos[[p]], na.rm=T, plot=F)
  S<-nrow(ll[[1]]$rsquare)
  res[p]<-(sum(ll[[1]]$rsquare, na.rm=T)*2)/(S*(S-1))

 } # fine for in p

 res

}


fst<-function(haplos, pop=c(), pairwise=TRUE, pop_ind=c()) {

 # calcola fst pairwise come descritto da estimatore Hudson,Slatkin,Maddison 1992 eq.3

 if(length(pop)>0) {
  nhaplos<-list()
  for (p in 1:length(pop)) {
   nhaplos[[p]]<-haplos[[pop[p]]]
  }
  haplos<-nhaplos
 }

 # quante subpop
 npop<-length(haplos)
 nsam<-rep(NA, npop)
 strhaplos<-list()
 for (p in 1:npop) {
  haplos[[p]]<-toupper(haplos[[p]])
  nsam[p]<-length(haplos[[p]])
  strhaplos[[p]]<-strsplit(haplos[[p]], split="")
 }

 # mean number of differences between different sequences from the same subpopulations
 ndiff<-list()
 Hw<-rep(NA, npop)
 for (p in 1:npop) {
  ndiff[[p]]<-matrix(NA, nrow=nsam[p], ncol=nsam[p])
  for (i in 1:(nsam[p]-1)) {
   for (j in (i+1):nsam[p]) {
#     hi<-strsplit(haplos[[p]][i],split="")[[1]]
#     hj<-strsplit(haplos[[p]][j],split="")[[1]]
    hi<-strhaplos[[p]][[i]]
    hj<-strhaplos[[p]][[j]]
#    ndiff[[p]][i,j]<-length(which(hi!=hj))
    ndiff[[p]][i,j]<-length(which(hi!=hj & hi!="N" & hj!="N"))
    Hw[p]<-mean(ndiff[[p]], na.rm=T)
   }
  }
 }

 # mean number of differences between sequences sampled from the two different subpopulations sampled
# npair<-list()
 Hb<-matrix(NA, nrow=npop, ncol=npop)
 for (p1 in 1:(npop-1)) {
  for (p2 in (p1+1):npop) {
   npair<-matrix(NA, nrow=nsam[p1], ncol=nsam[p2])
   for (i in 1:(nsam[p1])) {
    for (j in 1:nsam[p2]) {
#      hi<-strsplit(haplos[[p1]][i],split="")[[1]]
#      hj<-strsplit(haplos[[p2]][j],split="")[[1]]
     hi<-strhaplos[[p1]][[i]]
     hj<-strhaplos[[p2]][[j]]
#     npair[i,j]<-length(which(hi!=hj))
     npair[i,j]<-length(which(hi!=hj & hi!="N" & hj!="N"))
     Hb[p1,p2]<-mean(npair, na.rm=T)
    }
   }
  } # fine for p2
 } # fine for p1 

 # fst
 value<-matrix(NA, nrow=npop, ncol=npop)
 if (pairwise) {
  for (i in 1:(npop-1)) {
   for (j in (i+1):npop) {
 #   value[i,j]<-1-(mean(unlist(c(ndiff[[i]],ndiff[[j]])),na.rm=T)/Hb[i,j])
    value[i,j]<-1-(mean(c(Hw[[i]],Hw[[j]]))/Hb[i,j])
   }
  }
 # value[1,1]<-(1-(mean( unlist(ndiff), na.rm=T )/ mean( unlist(npair), na.rm=T)))
 }
 value[1,1]<-(1-(mean( Hw, na.rm=T )/ mean( Hb, na.rm=T)))
 colnames(value)<-rownames(value)<-seq(1,npop)

 perc<-perc5k<-c()

 # risultato
# if (length(pop_ind==0)) res_final<-value else res_final<-list(value, perc, perc5k)

 res_final<-list(value=value, perc=perc, perc5k=perc5k)

 res_final

} # fine function




















