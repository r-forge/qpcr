setGeneric("qPCRCorMatrix",
  function(qPCRSet)
  standardGeneric("qPCRCorMatrix")
)
setMethod("qPCRCorMatrix", signature = "qPCRSet", definition =
  function(qPCRSet)
  {
    matrix_maker <- vector()
    raw <- as.data.frame(exprs(qPCRSet))
    no_samples <- length(sampleNames(qPCRSet))
    for (i in sampleNames(qPCRSet))
      for (j in sampleNames(qPCRSet)) matrix_maker <- c(matrix_maker, cor(raw[, i], raw[, j], use="complete.obs"))
    correlation_matrix <- matrix(matrix_maker,no_samples, no_samples)
    colnames(correlation_matrix) <- gsub(" ", "", sampleNames(qPCRSet))
    row.names(correlation_matrix) <- gsub(" ", "", sampleNames(qPCRSet))
    correlation_matrix[upper.tri(correlation_matrix)] <- NA
    return(correlation_matrix)
  }
)
