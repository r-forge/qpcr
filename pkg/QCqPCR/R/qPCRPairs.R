setGeneric("qPCRPairs",
  function(qPCRBatch, writeToFile=FALSE, pairsToPlot="AllPairs")
  standardGeneric("qPCRPairs")
)
setMethod("qPCRPairs", signature = "qPCRBatch", definition =
  function(qPCRBatch, writeToFile, pairsToPlot) {
    if(pairsToPlot != "AllPlates") { 
        pairsToPlot <- combn(sampleNames(qPCRBatch),2)
    }
    else {
        pairsToPlot <- combn(pairsToPlot,2)
    }
    apply(pairsToPlot, 2, .plotPairs, qPCRBatch, writeToFile)
  }
)
.plotPairs <- function(samples, qPCRBatch, writeToFile) # plots graph between the 2 samples
{
  x <- samples[1]
  y <- samples[2]
  if (writeToFile) jpeg(filename = paste(x, "_", y, "_Pairs_Plot.jpeg", sep = ""))
  raw <- as.data.frame(exprs(qPCRBatch))
  plot(raw[, x], raw[, y], xlab = x, ylab = y, xlim = c(1, max(raw[, x], na.rm=TRUE)), ylim = c(1, max(raw[, y], na.rm=TRUE)))
  title(main <- c(x, "vs", y, "R^2 = ", cor(raw[, x] ,raw[, y], use <- "complete.obs")))
  abline(0, 1)
  if (writeToFile) dev.off()
}
