#setClass("normqPCRSet", # tSet is a matrix of all the detectors x different values for each hkg, with a vector of hkgs
#         representation = representation(
#           nSet = "data.frame", 
#           hkgs = "character"))

###################################
#    normalise qSet with 1 HKG    #
###################################

normaliseByHKG <- function(qSet, hkgs, design, verbose = FALSE){ # takes expression set, vector of housekeeping genes and a design matrix (like in limma)

    ##########
    # Use design 'matrix' to work out which is case and which is control
    logicase <- design == "case" 
    logicontrol <- design == "control"
    lenCase <- sum(logicase==TRUE)
    lenControl <- sum(logicontrol==TRUE)
    maxNAinCase <- lenCase - 2 # max number of NA values allowed for case
    maxNAinControl <- lenControl - 2 # max number of NA values allowed for control

    ##########

    normSet <- data.frame(as.data.frame(exprs(qSet))) # turn exprs component into a data frame
    tabFormat <- vector()
    laterColNames <- vector()

    for(hkg in hkgs) { # For each nominal housekeeping gene
        ##########
        # first see if it's suitable
        hkgCts <- as.numeric(exprs(qSet[hkg,])) # get the Ct values for given hkg
        hkgCtsCase <- hkgCts[logicase] # seperate for case and control
        hkgCtsControl <- hkgCts[logicontrol]
        # if unsuitable, stop the loop
        if(sum(is.na(hkgCtsCase)) > maxNAinCase) stop (hkg, " is unsuitable as a housekeeping gene because a value was only obtained for it once or less")
        if(sum(is.na(hkgCtsControl)) > maxNAinControl) stop (hkg, " is unsuitable as a housekeeping gene because a value was only obtained for it once or less")

        if(hkg %in% featureNames(qSet) == FALSE) stop (hkg," not found in file. Ensure entered housekeeping genes appear in the file")
        hkg <- gsub("-.+$","",hkg) # regexp to remove any rubbish from end of control gene spec
	laterColNames <- c(laterColNames, paste(hkg, "Control_mean", sep = "_"), paste(hkg, "Control_Sds", sep = "_"), paste(hkg, "Case_mean", sep = "_"), paste(hkg, "Case_Sds", sep = "_"), paste(hkg, "ddCt", sep = "_"), paste(hkg, "2^DDCt", sep = "_"), paste(hkg, "2^DDCt min", sep = "_"), paste(hkg, "2^DDCt max", sep = "_")) # get the column names for the tSet matrix
    }

    for(hkg in hkgs) { # for each housekeeping gene
        if(verbose) cat("For HKG: ", hkg, "\n")
        if(hkg %in% featureNames(qSet) == FALSE) stop (hkg," not found in file. Ensure entered housekeeping genes appear in the file")
        hkgCts <- as.numeric(exprs(qSet[hkg, ]))
        hkgCtsCase <- hkgCts[logicase]
        hkgCtsControl <- hkgCts[logicontrol]
        ##########
        # first see if it's suitable
        hkgCts <- as.numeric(exprs(qSet[hkg, ])) # get the Ct values for given hkg
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
        for (detector in featureNames(qSet)) {
            Cts <- as.numeric(exprs(qSet[detector,])) # the raw values for the detector
            Cts[Cts>38] <- NA
            CtsControl <- Cts[logicontrol]
            CtsCase <- Cts[logicase]
            if(sum(is.na(CtsControl)) > maxNAinControl) { # if we have under 2 values
                normControl <- NA
                v_dCtsMeanControl <- NA
                v_dCtsSdsControl <- NA
            }
            else {
                normControl <- mean(CtsControl - hkgCtsControl,na.rm=T) # otherwise work out stats
                v_dCtsMeanControl <- mean(CtsControl - hkgCtsControl,na.rm=T)
                v_dCtsSdsControl <- sd(CtsControl - hkgCtsControl,na.rm=T)
            }
            if(sum(is.na(CtsCase)) > maxNAinCase) { # if we have under 2 values
                normCase <- NA
                v_dCtsMeanCase <- NA
                v_dCtsSdsCase <- NA
            }
            else {
                normCase <- mean(CtsCase - hkgCtsCase,na.rm=T) # otherwise work out stats
                v_dCtsMeanCase <- mean(CtsCase - hkgCtsCase,na.rm=T)
                v_dCtsSdsCase <- sd(CtsCase - hkgCtsCase,na.rm=T)
            }
            deltaCtsMeanCase <- c(deltaCtsMeanCase,v_dCtsMeanCase)
            deltaCtsSdCase <- c(deltaCtsSdCase,v_dCtsSdsCase)
            deltaCtsMeanControl <- c(deltaCtsMeanControl,v_dCtsMeanControl)
            deltaCtsSdControl <- c(deltaCtsSdControl,v_dCtsSdsControl)

            if (is.na(normControl) & is.na(normCase)) {} else
            if (is.na(normControl)) {normControl = 38 - min(hkgCts)} else # will fall over if NAs in HKG
            if (is.na(normCase)) {normCase = 38 - min(hkgCts)}

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
    normSetColNames <- c(sampleNames(qSet),laterColNames) # make column names
    names(normSet) <- normSetColNames
    qSet@nSet <- normSet
    return(qSet)
}
