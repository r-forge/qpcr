setGeneric("combineTechReps",
  function(qPCRBatch)
  standardGeneric("combineTechReps")
)
setMethod("combineTechReps", signature = "qPCRBatch", definition =
  function(qPCRBatch) {
    expM <- exprs(qPCRBatch)
    origDetectors <- row.names(expM)
    if (FALSE %in% grepl("_TechReps", origDetectors)) stop("These are not tech reps")
    newDetectors <- unique(gsub("_TechReps.\\d","", origDetectors))
    NewExpM <- matrix(nrow = length(newDetectors), ncol = dim(expM)[2], dimnames = list(newDetectors,colnames(expM)))
    for(detector in newDetectors){
	cat(detector,"\t")
        cat(expM[gsub("_TechReps.\\d","", origDetectors) %in% origDetectors],"\t")
      dValues <- colMeans(expM[grepl(detector, origDetectors),],na.rm=TRUE)
	cat(dValues,"\n")
      NewExpM[detector,] <- dValues
    }
    NewExpM[is.na(NewExpM)] <- NA
    qPCRBatch <- new("qPCRBatch", exprs = NewExpM, phenoData = phenoData(qPCRBatch))
    return(qPCRBatch)
  }
)
