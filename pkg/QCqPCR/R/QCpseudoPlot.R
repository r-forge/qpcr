setGeneric("PseudoPlot",
  function(qPCRSet, plotType="Cts.Values")
  standardGeneric("PseudoPlot")
)

setMethod("PseudoPlot", signature = "qPCRSet", definition =
  function(qPCRSet, plotType) {
    ctsMat <- exprs(qPCRSet)
    ctsMat[is.na(ctsMat)] <- 40 # Cutoff point..Could change for different platforms
    orderMat <- exprs.well.order(qPCRSet)
    plateVec <- as.vector(gsub("-.*", "", orderMat))
    wellVec <- as.numeric(gsub(".*-", "", orderMat))

    if (plotType == "Cts.Values") {
      minVal <- round(min(ctsMat, na.rm=TRUE), 2)
      maxVal <- round(max(ctsMat, na.rm=TRUE), 2)

      for (plate in unique(plateVec)) { # for each plate
        plotTitle <- paste(plotType, "for plate:", plate)
        orderedVals <- ctsMat[plateVec == plate][order(wellVec[plateVec == plate])]
        plotMat <- matrix(orderedVals, nrow=16, byrow=TRUE)
        .plotCard(plotMat, plotTitle, minVal, maxVal)
        .wait()
      }
    }
    else if (plotType == "Plate.Residuals") {
      for (plate in unique(plateVec)) {
        plotTitle <- paste(plotType, "for plate:", plate)
        orderedVals <- ctsMat[plateVec == plate][order(wellVec[plateVec == plate])]
        meanSubbedVals <- abs(orderedVals - mean(orderedVals, na.rm=TRUE))
        plotMat <- matrix(meanSubbedVals, nrow=16, byrow=TRUE)
        minVal <- 0
        maxVal <- round(max(plotMat, na.rm=TRUE), 2)
        .plotCard(plotMat, plotTitle, minVal, maxVal)
        .wait()
      }
    }
    else if (plotType == "Detector.Residuals") {
      valMat <- abs(ctsMat - rowMeans(ctsMat, na.rm=TRUE)) # take the avg values from the Cts vals
      minVal <- 0
      maxVal <- round(max(valMat, na.rm=TRUE), 2)
      for (plate in unique(plateVec)) {
        title <- paste(plotType, "for plate:", plate)
        orderedVals <- valMat[plateVec == plate][order(wellVec[plateVec == plate])]
        plotMat <- matrix(orderedVals, nrow=16, byrow=TRUE)
        .plotCard(plotMat, title, minVal, maxVal)
        .wait()
      }
    }
    else if (plotType == "Well.Residuals") {
      averageWell <- vector(length = max(wellVec)) # Initialise a vector of the average Ct value
      for (well in 1:max(wellVec)) { # generate average well amounts
        wellChar <- as.character(well)
        averageWell[well] <- mean(ctsMat[wellVec == wellChar], na.rm=TRUE) # add the mean value for a given well
      }
      for (plate in unique(plateVec)) { # for each plate
        title <- paste(plotType, "for plate:", plate)
        orderedCts <- ctsMat[plateVec == plate][order(wellVec[plateVec == plate])]
        wellMeanSubbedCts <- abs(orderedCts - averageWell)
        plotMat <- matrix(wellMeanSubbedCts, nrow=16, byrow=TRUE)
        minVal <- round(min(plotMat, na.rm=TRUE), 2)
        maxVal <- round(max(plotMat, na.rm=TRUE), 2)
        .plotCard(plotMat, title, minVal, maxVal)
        .wait()
      }
    }
  }
)

.plotCard = function(plotMat, plotTitle, minVal, maxVal) 
{

  if(writeToFile == TRUE) {
    jpegTitle <- paste(plotTitle,".jpg",sep="")
    jpeg(jpegTitle)
  }

  background <- "black"
  n <- nrow(plotMat)
  m <- ncol(plotMat)
  rname <- 1:n
  cname <- 1:m
  rname <- as.character(rname)
  cname <- as.character(cname)
  plot.new()
  plot.window(c(0, m), c(0, n), asp = 1)
  xlabwidth <- max(strwidth(rname, cex = 1))
  ylabwidth <- max(strwidth(cname, cex = 1))
  ## set up an empty plot with the appropriate dimensions
  plot.window(c(-xlabwidth + 0.5, m + 0.5), c(0, n + 1 + ylabwidth), asp = 1, xlab="", ylab="")
  rect(0.5, 0.5, m + 0.5, n + 0.5, col = background)      ##background color
  title(plotTitle)
  ## Grid
  segments(rep(0.5, n + 1), 0.5 + 0:n, rep(m + 0.5, n + 1), 0.5 + 0:n, col = "gray")
  segments(0.5 + 0:m, rep(0.5, m + 1), 0.5 + 0:m, rep(n + 0.5, m), col = "gray")
  # now work out the colours to use
  rampSeq <- seq(minVal, maxVal, 0.01)
  col <- colorRampPalette(c("green","red"))(length(rampSeq))
  bg <- col[as.vector((plotMat - minVal) * 100 + 1)]
#   # and plot the symbols with colours
  symbols(rep(1:m, each = n), rep(n:1, m), add = TRUE, inches = F, circles = rep(0.4, (m*n)), bg = as.vector(bg))

  if(writeToFile == TRUE) {
    dev.off()
  }
}

.wait <- function() {
  par(ask=TRUE)
}