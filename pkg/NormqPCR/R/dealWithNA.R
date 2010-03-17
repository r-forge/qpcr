#######################################################
# replaceNAs - this replaces all NAs with a set number
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
###############################################################################################
# replaceAboveCutOff - this replaces anything above a given number with a (supplied) new value
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
################################################################################################################
# makeAllNAs - for each detector, if you have > a given number of NAs, then all values are all replaced with NA
# This means we can ignore any NAs in future calculations (since they can be dealt with using these functions
setGeneric("makeAllNAs",
  function(qPCRBatch, contrastM, sampleMaxM)
  standardGeneric("makeAllNAs")
)
setMethod("makeAllNAs", signature = "qPCRBatch", definition =
  function(qPCRBatch, contrastM, sampleMaxM) {
    expM <- exprs(qPCRBatch)
    for (detector in featureNames(qPCRBatch)) {
      for(phenotype in colnames(sampleMaxM)) {
        pColumns <- row.names(contrastM)[contrastM[,phenotype] == 1]
      if(sum(is.na(expM[detector,pColumns])) >= sampleMaxM[,phenotype]) expM[detector,pColumns] <- NA
      }
    }
    exprs(qPCRBatch) <- expM
    return(qPCRBatch)
  }
)
