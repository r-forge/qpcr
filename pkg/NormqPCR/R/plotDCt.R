plotDCt <- function(..., ddCtTable, detectors="", statCalc="arith") {
  if(detectors[1]!="") {
    ddCtTable <- ddCtTable[ddCtTable$ID %in% detectors, ,drop=FALSE]
  } else {
    ddCtTable <- ddCtTable
  }
  plotNames <- ddCtTable$ID
  plotCts <- ddCtTable[,c(2,4)]
  plotSds <- ddCtTable[,c(3,5)]
  plotCts <- sapply(plotCts, function(x) as.numeric(as.character(x)))
  plotSds <- sapply(plotSds, function(x) as.numeric(as.character(x)))

  if(statCalc == "arith") {
    plotU <- 2^(log2(plotCts) + plotSds)
    plotL <- 2^(log2(plotCts) - plotSds)
  } else {
    plotU <- plotCts  + plotSds
    plotL <- plotCts - plotSds
  }
#  print(plotCts)
  
#  cat(plotU,"\t")
#  cat(plotL,"\n")

  barplot2(..., height=t(plotCts), plot.ci=TRUE, ci.u=t(plotU), ci.l=t(plotL), beside=T)
}
