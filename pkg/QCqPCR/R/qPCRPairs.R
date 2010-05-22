setGeneric("qPCRPairs",
  function(qPCRBatch, plotType="Sample", writeToFile=FALSE, pairsToPlot="All")
  standardGeneric("qPCRPairs")
)
setMethod("qPCRPairs", signature = "qPCRBatch", definition =
  function(qPCRBatch, plotType, writeToFile, pairsToPlot) {
    if(plotType == "Sample") {
      if(pairsToPlot == "All") {
          pairsToPlot <- combn(sampleNames(qPCRBatch),2)
      }
      else {
cat("HEREIHOPE\n")
          pairsToPlot <- combn(pairsToPlot,2)
      }
      plotMat <- exprs(qPCRBatch)
      apply(pairsToPlot, 2, .plotPairs, plotMat, writeToFile)
    }
    else if (plotType == "Plate") {
      if(pairsToPlot == "All") {
          pairsToPlot <- combn(sampleNames(qPCRBatch),2)
      }
      else {
        plateVec <- as.vector(gsub("-.*", "", pairsToPlot)) # orderMat not defined!
        wellVec <- as.numeric(gsub(".*-", "", pairsToPlot))
        plotMat <- matrix(ncol = length(levels(as.factor(plateVec))), nrow = max(wellVec))
      }
        for(i in plateVec) {
          plotMat[,count] <-  exprs(qPCRBatch)[plateVec == i][wellVec[plateVec == i]]
          pairsToPlot <- combn(pairsToPlot,2)
          count <- count + 1
        }
 
      apply(pairsToPlot, 2, .plotPairs, plotMat, writeToFile)
    }
    else stop("incorrect plotType argument given")
}
)
.plotPairs <- function(samples, plotMat, writeToFile) # plots graph between the 2 samples
{
  x <- samples[1]
  y <- samples[2]
  if (writeToFile) jpeg(filename = paste(x, "_", y, "_Pairs_Plot.jpeg", sep = ""))
#  plotMat <- as.data.frame(exprs(qPCRBatch))
  plot(plotMat[, x], plotMat[, y], xlab = x, ylab = y, xlim = c(1, max(plotMat[, x], na.rm=TRUE)), ylim = c(1, max(plotMat[, y], na.rm=TRUE)))
  title(main = paste(x, "vs", y, "R^2 = ", cor(plotMat[, x] ,plotMat[, y], use = "complete.obs")))
  abline(0, 1)
  if (writeToFile) dev.off()
}
