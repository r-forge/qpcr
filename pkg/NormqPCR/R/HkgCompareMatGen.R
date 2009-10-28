#           nSet = "data.frame", 
#           hkgs = "character"))

###################################
#    normalise qPCRSet with 1 HKG    #
###################################

produceHKGSMat <- function(qPCRSet, hkgs, design, cutoff = 38, verbose = FALSE){ # takes qPCRSet, vector of housekeeping genes and a design vector to say what comparisons to make
    hkgs <- make.names(hkgs)
    ##########
    # Use design 'matrix' to work out which is case and which is control
    logicase <- design == "case" 
    logicontrol <- design == "control"
    lenCase <- sum(logicase==TRUE)
    lenControl <- sum(logicontrol==TRUE)
    maxNAinCaseHkg <- ceiling(lenCase/2) # max number of NA values allowed for case
    maxNAinControlHkg <- ceiling(lenControl/2) # max number of NA values allowed for control
    DiscardNACase <- floor(lenCase/4) # less or equal to this number of NAs we discard
    DiscardNAControl <- floor(lenControl/4) # max number of NA values allowed for control
    ##########

    normSet <- data.frame(as.data.frame(exprs(qPCRSet))) # turn exprs component into a data frame
    normSet[normSet > cutoff] <- NA
    tabFormat <- vector()
    laterColNames <- vector()

    for(hkg in hkgs) { # For each nominal housekeeping gene
        ##########
        # first see if it's suitable
        hkgCts <- as.numeric(exprs(qPCRSet[hkg,])) # get the Ct values for given hkg
        hkgCtsCase <- hkgCts[logicase] # seperate for case and control
        hkgCtsControl <- hkgCts[logicontrol]
        # if unsuitable, stop the loop
        if(sum(is.na(hkgCtsCase)) > maxNAinCaseHkg) stop (hkg, " is unsuitable as a housekeeping gene because a value was obtained for it less than", maxNAinCaseHkg, "times out of", lenControl, ".")
        if(sum(is.na(hkgCtsControl)) > maxNAinControlHkg) stop (hkg, " is unsuitable as a housekeeping gene because a value was obtained for it less than", maxNAinControlHkg, "times out of ", lenControl, ".")

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
        if(sum(is.na(hkgCtsCase)) > maxNAinCaseHkg) stop (hkg, " is unsuitable as a housekeeping gene because a value was only obtained for it once or less")
        if(sum(is.na(hkgCtsControl)) > maxNAinControlHkg) stop (hkg, " is unsuitable as a housekeeping gene because a value was only obtained for it once or less")
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
            Cts <- as.numeric(exprs(qPCRSet[detector,])) # the raw values for the detector
            Cts[Cts > cutoff] <- NA
            CtsControl <- Cts[logicontrol]
            CtsCase <- Cts[logicase] 
            if(sum(is.na(CtsControl)) > DiscardNAControl) { # if up to a quarter NAs
                normControl <- mean(CtsControl - hkgCtsControl, na.rm=T) # work out stats using existing values
                v_dCtsMeanControl <- mean(CtsControl - hkgCtsControl, na.rm=T)
                v_dCtsSdsControl <- sd(CtsControl - hkgCtsControl, na.rm=T)
            }
            else { # otherwise make all the NAs cutoff value and continue
                CtsControl[CtsControl == NA] <- cutoff
                hkgCtsControl[hkgCtsControl == NA] <- cutoff	
                normControl <- mean(CtsControl - hkgCtsControl)
                v_dCtsMeanControl <- mean(CtsControl - hkgCtsControl)
                v_dCtsSdsControl <- sd(CtsControl - hkgCtsControl)
            }
            if(sum(is.na(CtsCase)) > DiscardNACase) { # if up to a quarter NAs
                normCase <- mean(CtsCase - hkgCtsCase,na.rm=T) # work out stats using existing values
                v_dCtsMeanCase <- mean(CtsCase - hkgCtsCase,na.rm=T)
                v_dCtsSdsCase <- sd(CtsCase - hkgCtsCase,na.rm=T)
            }
            else {
                CtsCase[CtsControl == NA] <- cutoff
                hkgCtsCase[hkgCtsControl == NA] <- cutoff
                normCase <- mean(CtsCase - hkgCtsCase)
                v_dCtsMeanCase <- mean(CtsCase - hkgCtsCase)
                v_dCtsSdsCase <- sd(CtsCase - hkgCtsCase)
            }
            ######## Now add these values to the vectors
            deltaCtsMeanCase <- c(deltaCtsMeanCase,v_dCtsMeanCase)
            deltaCtsSdCase <- c(deltaCtsSdCase,v_dCtsSdsCase)
            deltaCtsMeanControl <- c(deltaCtsMeanControl,v_dCtsMeanControl)
            deltaCtsSdControl <- c(deltaCtsSdControl,v_dCtsSdsControl)

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
            deltadeltaCt <- c(deltadeltaCt,ddct) # add to vector
        }
	normSet <- data.frame(normSet,deltaCtsMeanControl,deltaCtsSdControl,deltaCtsMeanCase,deltaCtsSdCase,deltadeltaCt,twoDDCt,minBound,maxBound) # add the vectors together
    }
    normSetColNames <- c(sampleNames(qPCRSet),laterColNames) # make column names
    names(normSet) <- normSetColNames
    return(normSet)
}
