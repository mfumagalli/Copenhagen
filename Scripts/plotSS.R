
#

win=5e4

fout="Results/CHB.ss.pdf"

# fins=c("Results/TSI.thetas.pestPG", "Results/CHB.thetas.pestPG", "Results/PEL.thetas.pestPG")
fins="Results/CHB.thetas.pestPG"

theta1 <- read.table(fins[1], fill=TRUE)
#theta2 <- read.table(fins[2], fill=TRUE)
#theta3 <- read.table(fins[3], fill=TRUE)

xl=c(109e6, 110e6)

pdf(file=fout)

par(mfrow=c(2,2))

plot(x=theta1$V3, y=theta1$V4, pch=16, ty="l", col="blue", xlab="Chromosome 2", ylab="Theta (Watt.)", xlim=xl)
#points(x=theta2$V3, y=theta2$V4, pch=16, col="black", ty="l")
#points(x=theta3$V3, y=theta3$V4, pch=16, col="green", ty="l")
#legend("topleft", legend=c("TSI", "CHB", "PEL"), col=c("blue","black","green"), pch=16)

plot(x=theta1$V3, y=theta1$V5, pch=16, col="blue", xlab="Chromosome 2", ylab="Theta (Pi)", xlim=xl, ty="l")
#points(x=theta2$V3, y=theta2$V5, pch=16, col="black", ty="l")
#points(x=theta3$V3, y=theta3$V5, pch=16, col="green", ty="l")
#legend("topleft", legend=c("TSI", "CHB", "PEL"), col=c("blue","black","green"), pch=16)

plot(x=theta1$V3, y=theta1$V5-theta1$V4, pch=16, col="blue", xlab="Chromosome 2", ylab="Tajima's D (~Pi-Watt)", xlim=xl, ty="l")
#points(x=theta2$V3, y=theta2$V5-theta2$V4, pch=16, col="black", ty="l")
#points(x=theta3$V3, y=theta3$V5-theta3$V4, pch=16, col="green", ty="l")
#legend("topleft", legend=c("TSI", "CHB", "PEL"), col=c("blue","black","green"), pch=16)

plot(x=theta1$V3, y=theta1$V14, pch=16, col="blue", xlab="Chromosome 2", ylab="Nr. of sites", xlim=xl, ty="l")
#points(x=theta2$V3, y=theta2$V14, pch=16, col="black", ty="l")
#points(x=theta3$V3, y=theta3$V14, pch=16, col="green", ty="l")
#legend("bottomright", legend=c("TSI", "CHB", "PEL"), col=c("blue","black","green"), pch=16)

dev.off()


cat("Output file:", fout, "\n")




