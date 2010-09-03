plotDCt <- function(ddCtTable, detectors="", logFC = FALSE) {
  if(detectors[1]!="") {
    ddCtTable <- ddCtTable[ddCtTable$ID %in% detectors, ,drop=FALSE]
  } else {
    ddCtTable <- ddCtTable
  }
  plotNames <- ddCtTable$ID
  plotCts <- ddCtTable[,c("case","control")]
  plotSds <- ddCtTable[,c("case.sd","control.sd")]
  plotCts <- sapply(plotCts, function(x) as.numeric(as.character(x)))
  plotSds <- sapply(plotSds, function(x) as.numeric(as.character(x)))

#  plotCts <- plotTable[,c("case","control")]
#  plotSds <- plotTable[,c("case.sd","control.sd")]
  plotU <- plotCts + plotSds
#  plotL <- plotCts - plotTable[,c("case.sd","control.sd")]
#  plotU <- plotCts + plotSds
  plotL <- plotCts - plotSds
print(plotCts)

#cat(plotU,"\n")
#cat(plotL,"\n")
#  plotCts <- sapply(ddCtTable, function(x) as.numeric(as.character(x)))
#  plotCtsd
#  plotCtMax
#  plotCtMin

  barplot2(height=t(plotCts), plot.ci=TRUE, ci.u=t(plotU), ci.l=t(plotL), beside=T)
}
