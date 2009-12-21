setGeneric("deltaDeltaAvgCt",
  function(qPCRBatch, maxNACase=0, maxNAControl=0, hkg, contrastM, case, control)
  standardGeneric("deltaDeltaAvgCt")
)
setMethod("deltaDeltaAvgCt", signature = "qPCRBatch", definition =
  function(qPCRBatch, maxNACase, maxNAControl, hkg, contrastM, case, control) {
    hkg <- make.names(hkg)
    if(! hkg %in% featureNames(qPCRBatch)) stop("invalid housekeeping gene")
    if(sum(is.na(hkg)) > 0) warning(hkg, " May be a bad housekeeping gene to normalise with since it did not produce a reading ", sum(is.na(hkg)), "times out of", length(hkg), ".")

    case <- row.names(contrastM)[contrastM[,case] == 1]
    control <- row.names(contrastM)[contrastM[,control] == 1]
    expM <- exprs(qPCRBatch)
    caseM <- expM[,case]
    controlM <- expM[,control]
    hkgVCase <- caseM[hkg, ]
    hkgVControl <- controlM[hkg, ]

    if(length(hkgVCase) == 1 || length(hkgVControl) == 1) {
      meanHkgCase <- hkgVCase
      meanHkgControl <- hkgVControl
      sdHkgCase <- NA
      sdHkgControl <- NA
    }
    else {
      meanHkgCase <- geomMean(hkgVCase, na.rm=TRUE)
      meanHkgControl <- geomMean(hkgVControl, na.rm=TRUE)
      sdHkgCase <- sd(hkgVCase, na.rm=TRUE)
      sdHkgControl <- sd(hkgVControl, na.rm=TRUE)
    }
    if(is.na(meanHkgCase) || is.na(meanHkgControl)) stop("Need at least 1 non NA for the housekeeper")

    ddCts <- vector(length=length(featureNames(qPCRBatch)))
    minddCts <- vector(length=length(ddCts))
    maxddCts <- vector(length=length(ddCts))

    i <- 1
    for (detector in featureNames(qPCRBatch)) {
          VCase <- caseM[detector,]
          VControl <- controlM[detector,]
          meanCase <- geomMean(VCase, na.rm=TRUE)
          meanControl <- geomMean(VControl, na.rm=TRUE)
          sdCase <- sd(VCase, na.rm=TRUE)
          sdControl <- sd(VControl, na.rm=TRUE)

        if(sum(is.na(VCase)) > maxNACase) {
          meanCase <- NA
          sdCase <- NA
        }
        else {
          meanCase <- geomMean(VCase, na.rm=TRUE)
          sdCase <- sd(VCase, na.rm=TRUE)
        }
        if(sum(is.na(VControl)) > maxNAControl) {
          meanControl <- NA
          sdControl <- NA
        }
        else {
          meanControl <- geomMean(VControl, na.rm=TRUE)
          sdControl <- sd(VControl, na.rm=TRUE)
        }
        dCtCase <- meanCase - meanHkgCase
        dCtControl <- meanControl - meanHkgControl
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
          sdCt <- sqrt((sdControl^2) + (sdCase^2))
          minddCts[i] <- 2 ^ -(ddCt + sdCt)
          maxddCts[i] <- 2 ^ -(ddCt - sdCt)
          ddCts[i] <- 2^-ddCt
        }
        i <- i+1
    }
    return(cbind(featureNames(qPCRBatch),ddCts,minddCts,maxddCts))
  }
)
