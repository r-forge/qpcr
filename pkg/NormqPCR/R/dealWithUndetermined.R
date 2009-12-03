setGeneric("replaceNAs",
  function(qPCRBatch, newNA=40)
  standardGeneric("replaceNAs")
)
setMethod("replaceNAs", signature = "qPCRBatch", definition =
  function(qPCRBatch, newNA) {
    exprs(qPCRBatch)[is.na(exprs(qPCRBatch))] <- newNA
    return(qPCRBatch)
  }
)

setGeneric("replaceAboveCutOff",
  function(qPCRBatch, newVal=NA, cutOff=38)
  standardGeneric("replaceAboveCutOff")
)
setMethod("replaceAboveCutOff", signature = "qPCRBatch", definition =
  function(qPCRBatch, newVal, cutOff) {
    exprs(qPCRBatch)[exprs(qPCRBatch) > cutOff] <- newVal
    return(qPCRBatch)
  }
)
