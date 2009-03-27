#           nSet = "data.frame", 
#           hkgs = "character"))

###################################
#    normalise qPCRSet with 1 HKG    #
###################################

normaliseByHKG <- function(qPCRSet, hkgs, design, cutoff = 38, verbose = FALSE){ # takes qPCRSet, vector of housekeeping genes and a design matrix (like in limma)

    ##########
    # Use design 'matrix' to work out which is case and which is control
    logicase <- design == "case" 
    logicontrol <- design == "control"
    lenCase <- sum(logicase==TRUE)
    lenControl <- sum(logicontrol==TRUE)
    maxNAinCase <- ceiling(lenCase/2) # max number of NA values allowed for case
    maxNAinControl <- ceiling(lenControl/2) # max number of NA values allowed for control
cat("LENGTH OF CASE IS:", lenCase, ":\n")
cat("LENGTH OF CONTROL IS:", lenControl, ":\n")
cat("MAX NA IS :", maxNAinCase, ":\n")
cat("MAX NA IS :", maxNAinControl, ":\n")
    ##########

    normSet <- data.frame(as.data.frame(exprs(qPCRSet))) # turn exprs component into a data frame
    normSet[normSet > 38] <- NA
    tabFormat <- vector()
    laterColNames <- vector()

    for(hkg in hkgs) { # For each nominal housekeeping gene
        ##########
        # first see if it's suitable
        hkgCts <- as.numeric(exprs(qPCRSet[hkg,])) # get the Ct values for given hkg
        hkgCtsCase <- hkgCts[logicase] # seperate for case and control
        hkgCtsControl <- hkgCts[logicontrol]
        # if unsuitable, stop the loop
        if(sum(is.na(hkgCtsCase)) > maxNAinCase) stop (hkg, " is unsuitable as a housekeeping gene because a value was only obtained for it once or less")
        if(sum(is.na(hkgCtsControl)) > maxNAinControl) stop (hkg, " is unsuitable as a housekeeping gene because a value was only obtained for it once or less")

        if(hkg %in% featureNames(qPCRSet) == FALSE) stop (hkg," not found in file. Ensure entered housekeeping genes appear in the file")
        hkg <- gsub("-.+$","",hkg) # regexp to remove any rubbish from end of control gene spec
	laterColNames <- c(laterColNames, paste(hkg, "Control_mean", sep = "_"), paste(hkg, "Control_Sds", sep = "_"), paste(hkg, "Case_mean", sep = "_"), paste(hkg, "Case_Sds", sep = "_"), paste(hkg, "ddCt", sep = "_"), paste(hkg, "2^DDCt", sep = "_"), paste(hkg, "2^DDCt min", sep = "_"), paste(hkg, "2^DDCt max", sep = "_")) # get the column names for the tSet matrix
    }

    for(hkg in hkgs) { # for each housekeeping gene
        if(verbose) cat("For HKG: ", hkg, "\n")
        if(hkg %in% featureNames(qPCRSet) == FALSE) stop (hkg," not found in file. Ensure entered housekeeping genes appear in the file")
        hkgCts <- as.numeric(exprs(qPCRSet[hkg, ]))
        hkgCtsCase <- hkgCts[logicase]
        hkgCtsControl <- hkgCts[logicontrol]
        ##########
        # first see if it's suitable
        hkgCts <- as.numeric(exprs(qPCRSet[hkg, ])) # get the Ct values for given hkg
        hkgCtsCase <- hkgCts[logicase] # seperate for case and control
        hkgCtsControl <- hkgCts[logicontrol]
        # if unsuitable, stop the loop
        if(sum(is.na(hkgCtsCase)) > maxNAinCase) stop (hkg, " is unsuitable as a housekeeping gene because a value was only obtained for it once or less")
        if(sum(is.na(hkgCtsControl)) > maxNAinControl) stop (hkg, " is unsuitable as a housekeeping gene because a value was only obtained for it once or less")
        ##########
        # initialise vectors for different stats for the different detectors for this hkg and either case/control
        deltaCtsMeanCase <- vector()
        deltaCtsSdCase <- vector()
        deltaCtsMeanControl <- vector()
        deltaCtsSdControl <- vector()
        deltadeltaCt <- vector()
        twoDDCt <- vector()
        minBound <- vector()
        maxBound <- vector()
        ##########
        for (detector in featureNames(qPCRSet)) {
	cat("DetectorName :", detector, ":\n",file = "dbg.txt", append=T)
            Cts <- as.numeric(exprs(qPCRSet[detector,])) # the raw values for the detector
            Cts[Cts > cutoff] <- NA
            CtsControl <- Cts[logicontrol]
            CtsCase <- Cts[logicase]
	cat("number of NAs is :", sum(is.na(CtsControl)), ":\n", file="dbg.txt", append=T)
            if(sum(is.na(CtsControl)) > maxNAinControl) { # if we have under 2 values
                normControl <- NA
                v_dCtsMeanControl <- NA
                v_dCtsSdsControl <- NA
            }
            else {
                CtsControl[CtsControl == NA] <- 38
                hkgCtsControl[hkgCtsControl == NA] <- 38	
                normControl <- mean(CtsControl - hkgCtsControl,na.rm=T) # otherwise work out stats
                v_dCtsMeanControl <- mean(CtsControl - hkgCtsControl,na.rm=T)
                v_dCtsSdsControl <- sd(CtsControl - hkgCtsControl,na.rm=T)
            }
	cat("number of NAs is :", sum(is.na(CtsCase)), ":\n", file="dbg.txt", append=T)
            if(sum(is.na(CtsCase)) > maxNAinCase) { # if we have under 2 values
                normCase <- NA
                v_dCtsMeanCase <- NA
                v_dCtsSdsCase <- NA
            }
            else {
                CtsCase[CtsControl == NA] <- 38
                hkgCtsCase[hkgCtsControl == NA] <- 38
                normCase <- mean(CtsCase - hkgCtsCase,na.rm=T) # otherwise work out stats
                v_dCtsMeanCase <- mean(CtsCase - hkgCtsCase,na.rm=T)
                v_dCtsSdsCase <- sd(CtsCase - hkgCtsCase,na.rm=T)
            }
            deltaCtsMeanCase <- c(deltaCtsMeanCase,v_dCtsMeanCase)
            deltaCtsSdCase <- c(deltaCtsSdCase,v_dCtsSdsCase)
            deltaCtsMeanControl <- c(deltaCtsMeanControl,v_dCtsMeanControl)
            deltaCtsSdControl <- c(deltaCtsSdControl,v_dCtsSdsControl)

 #           if (is.na(normControl) & is.na(normCase)) {} else
 #           if (is.na(normControl)) {normControl = 38 - min(hkgCts)} else # will fall over if NAs in HKG
 #           if (is.na(normCase)) {normCase = 38 - min(hkgCts)}

            ddct <- normCase - normControl # log ratio

            twoDDCt <- c(twoDDCt,2^-ddct) # fold change

            if(typeof(v_dCtsSdsControl) == "double") {
                maxBound <- c(maxBound,2^-(ddct - v_dCtsSdsControl))
                minBound <- c(minBound,2^-(ddct + v_dCtsSdsControl))
            }
            else {
                minBound <- c(minBound,NA)
                maxBound <- c(maxBound,NA)
            }
            deltadeltaCt <- c(deltadeltaCt,ddct)
        }
	normSet <- data.frame(normSet,deltaCtsMeanControl,deltaCtsSdControl,deltaCtsMeanCase,deltaCtsSdCase,deltadeltaCt,twoDDCt,minBound,maxBound) # add the vectors together
    }
    normSetColNames <- c(sampleNames(qPCRSet),laterColNames) # make column names
    names(normSet) <- normSetColNames
    qPCRSet@nSet <- normSet
    qPCRSet@hkgs <- hkgs
    return(qPCRSet)
}
