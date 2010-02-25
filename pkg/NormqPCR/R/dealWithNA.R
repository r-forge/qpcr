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
  function(qPCRBatch, maxNACase, maxNAControl, contrastM, case, control)
  standardGeneric("makeAllNAs")
)
setMethod("makeAllNAs", signature = "qPCRBatch", definition =
  function(qPCRBatch, maxNACase, maxNAControl, contrastM, case, control) {
    case <- row.names(contrastM)[contrastM[,case] == 1]
    control <- row.names(contrastM)[contrastM[,control] == 1]
    expM <- exprs(qPCRBatch)
    caseM <- expM[,case]
    controlM <- expM[,control]

    for (detector in featureNames(qPCRBatch)) {
      if(sum(is.na(expM[detector,case])) > maxNACase) caseM[detector,case] <- NA
      if(sum(is.na(expM[detector,control])) > maxNAControl) controlM[detector,case] <- NA
    }
    return(qPCRBatch)
  }
)
