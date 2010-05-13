setGeneric("gM_ddAvgCt",
  function(qPCRBatch, maxNACase=0, maxNAControl=0, hkgs, contrastM, case, control)
  standardGeneric("gM_ddAvgCt")
)
setMethod("gM_ddAvgCt", signature = "qPCRBatch", definition =
  function(qPCRBatch, maxNACase, maxNAControl, hkgs, contrastM, case, control) {
    if(length(hkgs)<=1) stop("More than 1 houskeeping gene required")
    hkgs <- make.names(hkgs)
    if(FALSE %in% (hkgs %in% featureNames(qPCRBatch))) stop("invalid housekeeping gene")



    case <- row.names(contrastM)[contrastM[,case] == 1]
    control <- row.names(contrastM)[contrastM[,control] == 1]
    expM <- exprs(qPCRBatch)
    caseM <- expM[,case]
    controlM <- expM[,control]
    hkgM <- expM[hkgs,]
#return(hkgM)
    if (TRUE %in% apply(apply(hkgM, 1, is.na),2,sum)>0)  {
        warning("NAs present in housekeeping genes. NAs will be excluded when combining housekeepers to make pseudogenes")
        if (0 %in% apply(! apply(hkgM, 1, is.na),2,sum)) stop("Need at least 1 non NA for each housekeeper")
    }
    hkgMCase <- caseM[hkgs, ]
    hkgMControl <- controlM[hkgs, ]

hkgVCase <- apply(hkgMCase, 2, geomMean, na.rm=TRUE)
hkgVControl <- apply(hkgMControl, 2, geomMean, na.rm=TRUE)
    ##############################################
    # Here we use the geoMean function to find the geometric mean of the chosen housekeepers


    if(! FALSE %in% is.na(hkgVCase) || ! FALSE %in% is.na(hkgVControl)) stop("Need at least 1 non NA for the housekeeper")
    gMhkgControl <- geomMean(hkgVControl, na.rm=TRUE)
    gMhkgCase <- geomMean(hkgVCase, na.rm=TRUE)

    ddCts <- vector(length=length(featureNames(qPCRBatch)))
    minddCts <- vector(length=length(ddCts))
    maxddCts <- vector(length=length(ddCts))

    i <- 1
    for (detector in featureNames(qPCRBatch)) {
        VCase <- caseM[detector,]
        VControl <- controlM[detector,]
        if(! FALSE %in% is.na(VCase)) { 
          warning("No Detector for Case")
          dCtCase <- rep(NA, length(VCase))
          sdCase <- NA
        }
#        if(length(VCase) == 1) {
#          warning("Only one Detector for Case")
#          dCtCase <- VCase
#          sdCase <- NA
#        }
        else {
          dCtCase <- geomMean(VCase, na.rm=T) - gMhkgCase
          sdCase <- sd(VCase, na.rm=TRUE)
        }
#        if(length(VControl) == 1) {
#          warning("Only one Detector for Control")
#          dCtControl <- VControl
#        }
        if(! FALSE %in% is.na(VControl)) {            
          warning("No Detector for Control")
          dCtControl <- rep(NA, length(VControl))
          sdControl <- NA
        }
        else {
          dCtControl <- geomMean(VControl, na.rm=TRUE) - gMhkgControl
          sdControl <- sd(VControl, na.rm=TRUE)
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
          if(is.na(sdCase) | is.na(sdControl)) {
            minddCts[i] <- NA
            maxddCts[i] <- NA
          }
          else {
            sdCt <- sqrt((sdCase^2) + (sdControl^2))
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
