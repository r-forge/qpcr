setGeneric("plotVsHkg",
  function(qPCRBatch, hkgs, writeToFile=FALSE)
    standardGeneric("plotVsHkg")
)
setMethod("plotVsHkg", signature = "qPCRBatch", definition =
  function(qPCRBatch, hkgs, writeToFile)
  {
    cts <- exprs(qPCRBatch) # this refers to the actual data
#    hkgs <-  # these are the housekeeping genes the data has been normalised to
    if(FALSE %in% hkgs %in% featureNames(qPCRBatch)) stop("one or more housekeeper not found in exprs matrix")
    plotFrame <- row.names(cts)
    for (hkg in hkgs) {
        cat("wiv da", hkg, "\n")
#        hkg <- gsub("-.+$", "", hkg) # cut off the stuff from the detector's name after the - 
#        hkg2ddct <- paste(hkg, "_2^DDCt", sep = "")
#        ddct <- cts[, hkg2ddct]
#        plotFrame <- data.frame(plotFrame, ddct)
         
    }
    plotFrame <- plotFrame[, -1] # take off first column of detector names
    for (i in seq(hkgs)) { # now order by each hkg and print
        plotFrame <- plotFrame[order(plotFrame[, i]), ] # order by 1st column
        if (writeToFile) jpeg(file = paste("OrderedBy", hkgs[i], ".jpg", sep = ""))
        matplot(plotFrame, type = "l", pch = seq(hkgs), lty = seq(hkgs), log = "y", main = paste("Ordered By ", hkgs[i], sep = ""),xlab=paste("rank order of ",hkgs[i]))
        legend(1, tail(seq(hkgs), n = 1), hkgs, lty = seq(hkgs), col = seq(hkgs))
        if (writeToFile) dev.off()
    }
    return(plotFrame)
  }
)
