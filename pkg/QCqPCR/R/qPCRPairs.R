setGeneric("qPCRPairs",
  function(qPCRSet, writeToFile=FALSE)
  standardGeneric("qPCRPairs")
)
setMethod("qPCRPairs", signature = "qPCRSet", definition =
  function(qPCRSet, writeToFile=FALSE) {
    for (sampleA in sampleNames(qPCRSet))
      for (sampleB in sampleNames(qPCRSet)) .plotPairs(sampleA, sampleB, qPCRSet, writeToFile)
  }
)
.plotPairs <- function(x, y, qPCRSet, writeToFile) # plots graph between the 2 samples
{
  if (writeToFile) jpeg(filename = paste(x, "_", y, "_Pairs_Plot.jpeg", sep = ""))
  raw <- as.data.frame(exprs(qPCRSet))
  plot(raw[, x], raw[, y], xlab = x, ylab = y, xlim = c(1, max(raw[, x], na.rm=TRUE)), ylim = c(1, max(raw[, y], na.rm=TRUE)))
  title(main <- c(x, "vs", y, "R^2 = ", cor(raw[, x] ,raw[, y], use <- "complete.obs")))
  abline(0, 1)
  if (writeToFile) dev.off()
}
