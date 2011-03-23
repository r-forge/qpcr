setMethod("qPCRCorMatrix", signature = "qPCRBatch", definition =
  function(qPCRBatch)
  {
    matrixMaker <- vector()
    raw <- as.data.frame(exprs(qPCRBatch))
    noSamples <- length(sampleNames(qPCRBatch))
    for (i in sampleNames(qPCRBatch))
      for (j in sampleNames(qPCRBatch)) matrixMaker <- c(matrixMaker, cor(raw[, i], raw[, j], use="complete.obs"))
    correlationMatrix <- matrix(matrixMaker,noSamples, noSamples)
    colnames(correlationMatrix) <- gsub(" ", "", sampleNames(qPCRBatch))
    row.names(correlationMatrix) <- gsub(" ", "", sampleNames(qPCRBatch))
    correlationMatrix[upper.tri(correlationMatrix)] <- NA
    return(correlationMatrix)
  }
)
