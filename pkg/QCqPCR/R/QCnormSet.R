plotVsHkg = function(nqSet,saveToFile=FALSE)
{
    nset <- nqSet@nSet # this refers to the actual data
    hkgs <- nqSet@hkgs # these are the housekeeping genes the data has been normalised to
    plotFrame <- row.names(nset)
    for (hkg in hkgs) {
        hkg2ddct <- paste(hkg, "_2^DDCt", sep = "")
        ddct <- nset[, hkg2ddct]
        plotFrame <- data.frame(plotFrame, ddct)
    }
    plotFrame <- plotFrame[, -1] # take off first column of detector names
    for (i in seq(hkgs)) { # now order by each hkg and print
        plotFrame <- plotFrame[order(plotFrame[, i]),] # order by 1st column
        if (saveToFile) jpeg(file = paste("OrderedBy", hkgs[i], ".jpg", sep = ""))
        matplot(plotFrame, type = "l", pch = seq(hkgs), lty = seq(hkgs), log = "y", main = paste("Ordered By ", hkgs[i], sep = ""),xlab=paste("rank order of ",hkgs[i]))
        legend(1, tail(seq(hkgs), n = 1), hkgs, lty = seq(hkgs), col = seq(hkgs))
        if (saveToFile) dev.off()
    }
    return(plotFrame)
}
