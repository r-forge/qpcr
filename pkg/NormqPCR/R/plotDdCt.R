
plotDdCt <- function(..., ddCtTable, detectors="", logFC=TRUE) {
  if(detectors[1]!="") {
    ddCtTable <- ddCtTable[ddCtTable$ID %in% detectors, ,drop=FALSE]
  } else {
    ddCtTable <- ddCtTable
  }
  plotNames <- ddCtTable$ID
ds <- gsub("\\.\\w*$","",plotNames) 
  plotdDCt <- as.numeric(as.vector(ddCtTable[,6]))
  plotMin <- as.numeric(as.vector(ddCtTable[,7]))
  plotMax <- as.numeric(as.vector(ddCtTable[,8]))

plotdDCt[grep("\\+",as.vector(ddCtTable[,6]))] <- max(plotMax, na.rm=T)
plotdDCt[grep("\\-",as.vector(ddCtTable[,6]))] <- min(plotMin, na.rm=T)
bar.colours <- rep("blue",length(plotdDCt))
bar.colours[grep("\\+",as.vector(ddCtTable[,6]))] <- "red"
bar.colours[grep("\\-",as.vector(ddCtTable[,6]))] <- "red"


  if(logFC == TRUE) {
    plotdDCt <- log2(plotdDCt)
    plotMin <- log2(plotMin)
    plotMax <- log2(plotMax)
  }
  barplot2(..., names=ds, height=plotdDCt, plot.ci=TRUE, ci.u=plotMax, ci.l=plotMin, las=2, col=bar.colours, ylab="Log2 Fold Change")
}


