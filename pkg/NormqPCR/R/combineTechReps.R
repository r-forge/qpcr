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
    for(detector in newDetectors){ # for each detector
      dValues <- colMeans(expM[gsub("_TechReps.\\d","", origDetectors) %in% detector,],na.rm=TRUE) # find everything that begins with new detector, and find the mean
      NewExpM[detector,] <- dValues # add values to the matrix
    }
    NewExpM[is.na(NewExpM)] <- NA # make NAs real NAs
    qPCRBatch <- new("qPCRBatch", exprs = NewExpM, phenoData = phenoData(qPCRBatch))
    return(qPCRBatch)
  }
)
