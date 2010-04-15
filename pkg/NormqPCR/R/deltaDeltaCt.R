setGeneric("deltaDeltaCt",
  function(qPCRBatch, maxNACase=0, maxNAControl=0, hkg, contrastM, case, control)
  standardGeneric("deltaDeltaCt")
)
setMethod("deltaDeltaCt", signature = "qPCRBatch", definition =
  function(qPCRBatch, maxNACase, maxNAControl, hkg, contrastM, case, control) {
    hkg <- make.names(hkg)
    if(! hkg %in% featureNames(qPCRBatch)) stop("invalid housekeeping gene")
    if(sum(is.na(hkg)) > 0) warning(hkg, " May be a bad housekeeping gene to normalise with since it did not produce a reading ", sum(is.na(hkg)), "times out of", length(hkg), ". deltaDeltaAvgCt() function might give more robust results")

    case <- row.names(contrastM)[contrastM[,case] == 1]
    control <- row.names(contrastM)[contrastM[,control] == 1]
    expM <- exprs(qPCRBatch)
    caseM <- expM[,case]
    controlM <- expM[,control]
    hkgVCase <- caseM[hkg, ]
    hkgVControl <- controlM[hkg, ]

    if(! FALSE %in% is.na(hkgVCase) || ! FALSE %in% is.na(hkgVControl)) stop("Need at least 1 non NA for the housekeeper")

    ddCts <- vector(length=length(featureNames(qPCRBatch)))
    minddCts <- vector(length=length(ddCts))
    maxddCts <- vector(length=length(ddCts))

    i <- 1
    for (detector in featureNames(qPCRBatch)) {
        VCase <- caseM[detector,]
        VControl <- controlM[detector,]
#stop("length of case is:",VCase,"_",length(VCase))
        if(length(VCase) == 1) {
          warning("Only one Detector for Control")
          dCtCase <- VCase
          sdCase <- NA
        } else if(! FALSE %in% is.na(VCase)) {
          warning("No Detector for Case")
          dCtCase <- rep(NA, length = VCase)
          dCtControl <- NA
        } else {
          dCtCase <- geomMean(VCase, na.rm=TRUE) - geomMean(hkgVCase, na.rm=TRUE)
          sdCase <- sd(VCase - hkgVCase, na.rm=TRUE)
        }

        if(length(VControl) == 1) {
          warning("Only one Detector for Control")
          dCtControl <- VControl
          sdControl <- NA
        } else if(! FALSE %in% is.na(VControl)) {
          warning("No Detector for Control")
          dCtControl <- rep(NA, length = VControl)
          sdControl <- NA
        } else {
          dCtControl <- geomMean(VControl, na.rm=TRUE) - geomMean(hkgVControl, na.rm=TRUE)
        }
        if(sum(is.na(VCase)) > maxNACase) {
          dCtCase <- NA
        }
        if(sum(is.na(VControl)) > maxNAControl) {
          dCtControl <- NA
        }
        ddCt <- (dCtCase - dCtControl)
        
        if(is.na(ddCt)) {
          if(is.na(dCtCase) && ! is.na(dCtControl)) ddCt <- "-"
          else if(is.na(dCtControl) && ! is.na(dCtCase)) ddCt <- "+"
          else if(is.na(dCtControl) && is.na(dCtCase)) ddCt <- NA
          minddCts[i] <- NA
          maxddCts[i] <- NA
          ddCts[i] <- ddCt
        }
        else {
          if(is.na(sdCase)) {
            minddCts[i] <- NA
            maxddCts[i] <- NA
          }
          else {
            sdCt <- sqrt(sdCase)
            minddCts[i] <- 2 ^ -(ddCt + sdCt)
            maxddCts[i] <- 2 ^ -(ddCt - sdCt)
            ddCts[i] <- 2^-ddCt
          }
        }
        i <- i+1
    }
    return(cbind(featureNames(qPCRBatch),ddCts,minddCts,maxddCts))
  }
)
