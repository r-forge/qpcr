#library(gplots)
plotDdCt <- function(qPCRSet, ddCtTable, detectors="", logFC = FALSE) {
#	no.plots <- length(detectors)
cat("start")
#cat(eval(samples))

#expM

if(detectors!="") {
#cat("WEN")
    bbb <- ddCtTable[ddCtTable$ID %in% detectors,]
    qPCRSet <- qPCRSet[detectors]
#    expM <- matrix(exprs(qPCRSet)[detectors,], ncol=8)
#    rownames(expM) <- detectors
    expM <- exprs(qPCRSet)
} else {
cat("try")
    bbb <- ddCtTable
    expM <- exprs(qPCRSet)
}
cat("out")
#print(bbb)
cat("YEAH")
ddCts <- as.vector(bbb$ddCt)
ddCtsNames <- as.vector(bbb$ID)
ddCtsMax <- log2(as.numeric(as.vector(bbb$ddCt.max)))
ddCtsMin <- log2(as.numeric(as.vector(bbb$ddCt.min)))

maxo <- ceiling(max(ddCtsMax,na.rm=T))
mino <- floor(min(ddCtsMin,na.rm=T))
nonNAs.back <- log2(as.numeric(ddCts))
nonNAs.back[ddCts == "+"] <- maxo
nonNAs.back[ddCts == "-"] <- mino
ddCts <- nonNAs.back
op <- par(las=2)
bp <- barplot2(ddCts, plot.ci=TRUE, ci.l=ddCtsMin, ci.u=ddCtsMax, names.arg=ddCtsNames, ylab=expression("Fold change relative to control ("~Delta~Delta~"Cti)"))
par(op)
#for(i in rownames(qq)) { t<-c(t,(ww[i,]));points(rep(j[m],4),ww[i,]);m=m+1 }   
#qq <- qPCRSet
#cat("bbbbbbl")
#return(expM)
expM <- log2(expM)[,c(3:4,7:8)] - log2(expM)[,c(1:2,5:6)]
if(length(detectors) == 1) expM <- matrix(expM, ncol=4)
#cat(featureNames(qPCRSet))
#rownames(expM) <- featureNames(qPCRSet)
#print(expM)
#return(expM)
source("/home/bsm/jperkins/qpcr/pkg/NormqPCR/R/plotDdCt.R")
j<-1; for(i in 1:length(featureNames(qPCRSet))) {cat(i);cat(bp[j]);print(expM[i,]); points(rep(bp[j],4),expM[i,]);j=j+1 }

}
