setGeneric("qPCRBoxWhisker",
  function(qPCRSet, cutOff = 40, statType = "parametric", plotType = "sample")
    standardGeneric("qPCRBoxWhisker")
)
setMethod("qPCRBoxWhisker", signature = "qPCRSet", definition =
  function(qPCRSet, cutOff, statType, plotType) {

    if (statType == "parametric" || statType == "non-parametric") {
    }
    else {
        stop("Invalid statType argument, please use \"parametric\" or \"non-parametric\"")
    }
    ctsMat <- exprs(qPCRSet)
    ctsMat[is.na(ctsMat)] <- cutOff # Cutoff point..Could change for different platforms
    ctsMat[ctsMat > cutOff] <- cutOff
    orderMat <- exprs.well.order(qPCRSet)
    plateVec <- as.vector(gsub("-.*", "", orderMat))
    wellVec <- as.numeric(gsub(".*-", "", orderMat))

    if(plotType == "sample") {
        valueVector <- vector(length = length(ctsMat))
        sampleNameVec <- vector(length = length(ctsMat))
#BWVector <- vector(length = 100000)
        i <- 1
        sampleLength <- dim(ctsMat)[1]
cat("LENGTH IS:",sampleLength,"\n")
        for(sample in colnames(ctsMat)) {
cat("SAMPLE IS",sample,"\n")
            cat("I IS", i, "\n")
            valueVector[i:(i + sampleLength - 1)] <- ctsMat[,colnames(ctsMat) %in% sample]
            sampleNameVec[i:(i + sampleLength -1)] <- rep(sample, sampleLength)
            i <- i + sampleLength
        }
        cat(sampleNameVec,"\n\n")
#return(sampleNameVec)
        sampleNameVec <- factor(sampleNameVec)
#        cat(sampleNameVec,"\n\n")
        cat(valueVector,"\n")
        graphMaker <- cbind(valueVector, sampleNameVec)
        boxplot(valueVector~sampleNameVec, data = graphMaker, names=colnames(ctsMat), col="bisque")
    }
#    if(plotType == "well")
#    if(plotType == "plate")
    }
)
