setGeneric("plotVsHkg",
  function(qPCRBatch, hkgs, writeToFile=FALSE)
    standardGeneric("plotVsHkg")
)


setMethod("plotVsHkg", signature = "qPCRBatch", definition =
  function(qPCRBatch, hkgs, writeToFile)
  {
    cts <- exprs(qPCRBatch) # this refers to the actual data
    hkgs <- make.names(hkgs)

    if(FALSE %in% hkgs %in% featureNames(qPCRBatch)) stop("one or more housekeeper not found in exprs matrix")
    plotFrame <- matrix(nrow=length(featureNames(qPCRBatch)),ncol=length(hkgs), dimnames = list(featureNames(qPCRBatch), hkgs))
    for (hkg in hkgs) {
        cat("housekeeping genes", hkg, "\n\n")

	dCts <- deltaCt(qPCRBatch = qPCRBatch, hkgs = hkg, calc="arith")
        for(detector in featureNames(qPCRBatch)) {
#cat("\n\n\n")
#	print(dCts)
#	cat(exprs(qPCRBatch)[detector,],"\t")
	plotFrame[detector,hkg] <- mean(exprs(dCts)[detector,],na.rm=TRUE)
#cat(mean(exprs(qPCRBatch)[detector,],na.rm=TRUE),"\t")
#cat(plotFrame[detector,hkg],"\n")
#        hkg <- gsub("-.+$", "", hkg) # cut off the stuff from the detector's name after the - 
#        hkg2ddct <- paste(hkg, "_2^DDCt", sep = "")
#        ddct <- cts[, hkg2ddct]
#        plotFrame[,hkgs] <-  <- data.frame(plotFrame, ddct)
         }
#cat("\n")
    }

#    for (hkg in hkgs[1]) {
#        for(detectors in featureNames(qPCRBatch)) {
#	    cat(plotFrame[detector,hkg])

#        }
#    }
#stop()
#    plotFrame <- plotFrame[, -1] # take off first column of detector names
   for (hkg in hkgs) { # now order by each hkg and print
        ord.plotFrame <- plotFrame[order(plotFrame[, hkg]), ] # order by 1st column
        if (writeToFile) jpeg(file = paste("mean.deltaCt.ordered.by.", hkg, ".jpg", sep = ""))


        matplot(ord.plotFrame, type = "l", pch = seq(hkgs), lty = seq(hkgs), main = paste("Ordered By hkg ", hkg, sep = ""),xlab=paste("rank order of ",hkg))
        legend(ceiling(1/2*length(featureNames(qPCRBatch))), min(ord.plotFrame,na.rm=TRUE)+5, hkgs, lty = seq(hkgs), col = seq(hkgs))
        if (writeToFile) {
          dev.off()
        } else {
          .wait()
        }
    }
    return(plotFrame)
  }
)

.wait <- function() {
  par(ask=TRUE)
}
