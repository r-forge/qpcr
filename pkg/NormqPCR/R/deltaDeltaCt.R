setGeneric("deltaDeltaCt",
  function(qPCRBatch, maxNACase=0, maxNAControl=0, hkgs, contrastM, case, control, paired=TRUE, combineHkg=FALSE)
  standardGeneric("deltaDeltaCt")
)
setMethod("deltaDeltaCt", signature = "qPCRBatch", definition =
  function(qPCRBatch, maxNACase, maxNAControl, hkgs, contrastM, case, control, paired, combineHkg) {
    hkgs <- make.names(hkgs)
    if(combineHkg == TRUE) {
	if(length(hkgs) == 1) stop("Not enough hkgs given")
    }
#    for(hkg in hkgs) {
#    if(! hkg %in% featureNames(qPCRBatch)) stop("invalid housekeeping gene")
    if(FALSE %in% (hkgs %in% featureNames(qPCRBatch))) stop("invalid housekeeping gene given")
    for(hkg in hkgs){
        if(sum(is.na(hkg)) > 0) warning(hkg, " May be a bad housekeeping gene to normalise with since it did not produce a reading ", sum(is.na(hkg)), "times out of", length(hkg))
    }
    case <- row.names(contrastM)[contrastM[,case] == 1]
    control <- row.names(contrastM)[contrastM[,control] == 1]
    expM <- exprs(qPCRBatch)
    caseM <- expM[,case]
    controlM <- expM[,control]

#    hkgMCase <- caseM[hkgs, ]
#    hkgMControl <- controlM[hkgs, ]
    if(combineHkg == TRUE) {
	hkgMCase <- caseM[hkgs, ]
        hkgMControl <- controlM[hkgs, ]
	hkgVCase <- apply(hkgMCase, 2, geomMean, na.rm=TRUE)
	hkgVControl <- apply(hkgMControl, 2, geomMean, na.rm=TRUE)
    } else {
        hkg <- hkgs[1]
    }

    hkgVCase <- caseM[hkg, ]
    hkgVControl <- controlM[hkg, ]


    sdHkgCase <- sd(hkgVCase, na.rm=TRUE)
    sdHkgControl <- sd(hkgVControl, na.rm=TRUE)

    if(! FALSE %in% is.na(hkgVCase) || ! FALSE %in% is.na(hkgVControl)) stop("Need at least 1 non NA for the housekeeper")

    ddCts <- vector(length=length(featureNames(qPCRBatch)))
    dCtCases <- vector(length=length(featureNames(qPCRBatch)))
    dCtControls <- vector(length=length(featureNames(qPCRBatch)))
    sdCtCases <- vector(length=length(featureNames(qPCRBatch)))
    sdCtControls <- vector(length=length(featureNames(qPCRBatch)))
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
	  if (paired == TRUE) {
	    sdCase <- sd(VCase - hkgVCase, na.rm=TRUE)
	  } else  {
	    sdCase <- sqrt(sd(VCase, na.rm=TRUE)^2 + sdHkgCase^2)
	  }
        }

        if(length(VControl) == 1) {
          warning("Only one Detector for Control")
          dCtControl <- VControl
        } else if(! FALSE %in% is.na(VControl)) {
          warning("No Detector for Control")
          dCtControl <- rep(NA, length = VControl)
        } else {
          dCtControl <- geomMean(VControl, na.rm=TRUE) - geomMean(hkgVControl, na.rm=TRUE)
          if (paired == TRUE) {
            sdControl <- sd(VControl - hkgVControl, na.rm=TRUE)
          } else  {
            sdControl <- sqrt(sd(VControl, na.rm=TRUE)^2 + sdHkgControl^2)
          }
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
#            sdCt <- sqrt(sdCase)
            minddCts[i] <- 2 ^ -(ddCt + sdCase)
            maxddCts[i] <- 2 ^ -(ddCt - sdCase)
            ddCts[i] <- 2^-ddCt
          }
        }
	dCtCases[i] <- 2^-dCtCase
	sdCtCases[i] <- sdCase
	dCtControls[i] <- 2^-dCtControl
	sdCtControls[i] <- sdControl
        i <- i+1
    }
    return(cbind(featureNames(qPCRBatch),dCtCases,sdCtCases,dCtControls,sdCtControls,ddCts,minddCts,maxddCts))
  }
)
