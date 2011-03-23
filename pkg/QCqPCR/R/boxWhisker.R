setMethod("qPCRBoxWhisker", signature = "qPCRBatch", definition =
  function(qPCRBatch, cutOff = 40, statType = "parametric", plotType = "sample") {

    if (statType == "parametric" || statType == "non-parametric") {
    }
    else {
        stop("Invalid statType argument, please use \"parametric\" or \"non-parametric\"")
    }
    ctsMat <- exprs(qPCRBatch)
    ctsMat[is.na(ctsMat)] <- cutOff # Cutoff point..Could change for different platforms
    ctsMat[ctsMat > cutOff] <- cutOff
    orderMat <- exprs.well.order(qPCRBatch)
    plateVec <- as.vector(gsub("-.*", "", orderMat))
    wellVec <- as.numeric(gsub(".*-", "", orderMat))

    if(plotType == "sample") {
        valueVector <- vector(length = length(ctsMat))
        sampleNameVec <- vector(length = length(ctsMat))
        i <- 1
        sampleLength <- dim(ctsMat)[1]
cat("LENGTH IS:",sampleLength,"\n")
        for(sample in colnames(ctsMat)) {
cat("SAMPLE IS",sample,"\n")
            cat("I IS ", i, "\n")
            valueVector[i:(i + sampleLength - 1)] <- ctsMat[,colnames(ctsMat) %in% sample]
            sampleNameVec[i:(i + sampleLength -1)] <- rep(sample, sampleLength)
            i <- i + sampleLength
        }
        cat(sampleNameVec,"\n\n")
        sampleNameVec <- factor(sampleNameVec)
        cat(valueVector,"\n")
        graphMaker <- cbind(valueVector, sampleNameVec)
        boxplot(valueVector~sampleNameVec, data = graphMaker, names=colnames(ctsMat), col="bisque")
    }
    if(plotType == "plate") {
        valueVector <- vector(length = length(ctsMat))
        plateNameVec <- vector(length = length(ctsMat))
        i <- 1
        plateLength <- max(wellVec)

        cat("\nplateLength is:",plateLength,"\n")

        for(plate in levels(factor(plateVec))) {
            cat("I is ", i, "\n")
            valueVector[i:(i + plateLength - 1)] <- ctsMat[plateVec == plate]
            plateNameVec[i:(i + plateLength -1)] <- rep(plate, plateLength)
            i <- i + plateLength
        }
        plateNameVec <- factor(plateNameVec)
        graphMaker <- cbind(valueVector, plateNameVec)
        print(graphMaker)
        boxplot(valueVector~plateNameVec, data = graphMaker, names=levels(plateVec), col="bisque")
    }
  }
)
