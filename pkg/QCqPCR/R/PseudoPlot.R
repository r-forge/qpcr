setGeneric("PseudoPlot",
  function(qPCRBatch, plotType="Cts.Values", writeToFile=FALSE, cutOff=NA, statType="parametric", plateToPlot="AllPlates")
    standardGeneric("PseudoPlot")
)
setMethod("PseudoPlot", signature = "qPCRBatch", definition =
  function(qPCRBatch, plotType, writeToFile, cutOff, statType, plateToPlot) {
# CHECKING - IS THERE A CLEVERER WAY TO DO THIS?
    if (statType == "parametric"
      || statType == "non-parametric") {
    }
    else {
        stop("Invalid statType argument, please use \"parametric\" or \"non-parametric\"")
    }
    ctsMat <- exprs(qPCRBatch)
    if(! is.na(cutOff)) {
      ctsMat[is.na(ctsMat)] <- cutOff # Cutoff point..Could change for different platforms
      ctsMat[ctsMat > cutOff] <- cutOff
    }
    else {
      warning("No cutOff value given, if you are calculating residuals, the program it will crash out ungracefully")
    }
    orderMat <- exprs.well.order(qPCRBatch)
    plateVec <- as.vector(gsub("-.*", "", orderMat))
    whichPlates <- unique(plateVec)
    whichPlates <- sort(plateVec)
    if(plateToPlot != "AllPlates") whichPlates <- plateToPlot
    wellVec <- as.numeric(gsub(".*-", "", orderMat))

    if (plotType == "Cts.Values") {
      minVal <- 0
      maxVal <- round(max(ctsMat, na.rm=TRUE), 2)
      for (plate in whichPlates) { # for each plate
        plotTitle <- paste(plotType, "for plate:", plate)
        orderedVals <- ctsMat[plateVec == plate][order(wellVec[plateVec == plate])]
        plotMat <- matrix(orderedVals, nrow=16, byrow=TRUE)
        .plotCardRaw(plotMat, plotTitle, minVal, max(plotMat, na.rm=TRUE), writeToFile)
      }
    }
    else if (plotType == "Plate.Residuals") {
      for (plate in whichPlates) {
        plotTitle <- paste(plotType, "for plate:", plate)
        orderedVals <- ctsMat[plateVec == plate][order(wellVec[plateVec == plate])]
        if (statType == "parametric") {
          plateResidual <- sd(as.vector(orderedVals), na.rm = TRUE)
          plateTotalDispersion <-  orderedVals - mean(orderedVals, na.rm = TRUE)
        }
        else {
          plateResidual <- mad(as.vector(orderedVals), na.rm = TRUE)
          plateTotalDispersion <-  orderedVals - median(orderedVals, na.rm = TRUE)
        }
        plateResidualsMat = plateTotalDispersion / plateResidual
        plotMat <- matrix(plateResidualsMat, nrow=16, byrow=TRUE)
        .plotCardStats(plotMat, plotTitle, writeToFile, statType)
      }
    }
    else if (plotType == "Detector.Residuals") {
      if (statType == "parametric") {  
         totalMat <- ctsMat - rowMeans(ctsMat, na.rm=TRUE) # take the avg values from the Cts vals
         residVec <- apply(ctsMat, 1, sd, na.rm=TRUE) # take the sds for each row
      }
      else {
         totalMat <- ctsMat - rowMedians(ctsMat, na.rm=TRUE) # take the avg values from the Cts vals
         residVec <- apply(ctsMat, 1, mad, na.rm=TRUE) # take the mads for each row
      }
      valMat <-  totalMat /  residVec # now divide to get the results in terms of SDs/MADs from mean
      valMat[is.na(valMat)] <- 0 # bit cludgey - deals with when we have a 0 / 0 calculations
      for (plate in whichPlates) { # now we must order and plot the new values by plate
        plotTitle <- paste(plotType, "for plate:", plate)
        orderedVals <- valMat[plateVec == plate][order(wellVec[plateVec == plate])]
        plotMat <- matrix(orderedVals, nrow=16, byrow=TRUE)
        .plotCardStats(plotMat, plotTitle, writeToFile, statType)
      }
    }
    else if (plotType == "Well.Residuals") { # By well we mean position of the well on the card
      averageWell <- vector(length = max(wellVec)) # Initialise a vector of the average Ct value
      residWell <- vector(length = max(wellVec))
      for (well in 1:max(wellVec)) { # generate average well amounts and resids as a background - stops having to generate on the fly
        wellChar <- as.character(well)
        if (statType == "parametric") {
          averageWell[well] <- mean(ctsMat[wellVec == wellChar], na.rm=TRUE) # add the mean value for a given well
          residWell[well] <- sd(ctsMat[wellVec == wellChar], na.rm=TRUE) # add the SD value for a given well
        }
        else if (statType == "non-parametric") {
          averageWell[well] <- median(ctsMat[wellVec == wellChar], na.rm=TRUE) # add the mean value for a given well
          residWell[well] <- mad(ctsMat[wellVec == wellChar], na.rm=TRUE) # add the SD value for a given well
        }
      }
      for (plate in whichPlates) { # for each plate
        plotTitle <- paste(plotType, "for plate:", plate)
        orderedCts <- ctsMat[plateVec == plate][order(wellVec[plateVec == plate])]
        totalVec <- orderedCts - averageWell
        valMat <- totalVec / residWell
        valMat[is.na(valMat)] <- 0 # bit cludgey - deals with when we have a 0 / 0 calculations
        plotMat <- matrix(valMat, nrow=16, byrow=TRUE)
        .plotCardStats(plotMat, plotTitle, writeToFile, statType)
      }
    }
    else {
      stop("Invalid Plot Type")
    }
  }
)

.plotCardRaw = function(plotMat, plotTitle, minVal, maxVal, writeToFile)
{
  if(writeToFile == TRUE) {
    jpegTitle <- paste(plotTitle,".jpg",sep="")
    jpeg(jpegTitle)
  }
  myPal <- brewer.pal(5, "RdYlGn")
  myCol <- colorRampPalette(myPal[5:1])(128)
  background <- "black"
  n <- nrow(plotMat)
  m <- ncol(plotMat)
  rname <- 1:n
  cname <- 1:m
  rname <- as.character(rname)
  cname <- as.character(cname)

  layout(matrix(c(1, 2), 2, 1), heights = c(4,1))
  op <- par(mar = c(0, 4, 4, 4))

  plot.new()
  plot.window(c(0, m), c(0, n), asp = 1)
  xlabwidth <- max(strwidth(rname, cex = 1))
  ylabwidth <- max(strwidth(cname, cex = 1))
  plot.window(c(-xlabwidth + 0.5, m + 0.5), c(0, n + 1 + ylabwidth), asp = 1, xlab="", ylab="")
  rect(0.5, 0.5, m + 0.5, n + 0.5, col = background) #background color
  title(plotTitle)
  segments(rep(0.5, n + 1), 0.5 + 0:n, rep(m + 0.5, n + 1), 0.5 + 0:n, col = "gray")
  segments(0.5 + 0:m, rep(0.5, m + 1), 0.5 + 0:m, rep(n + 0.5, m), col = "gray")
  bg <- myCol[plotMat * (128/maxVal)]
  symbols(rep(1:m, each = n), rep(n:1, m), add = TRUE, inches = F, circles = rep(0.4, (m*n)), bg = as.vector(bg))
  x.bar <- seq(from = minVal, to = maxVal, length = 128)
  par(mar = c(5.1, 4.1, 1, 2))
  image(x.bar, 1, matrix(x.bar, length(x.bar), 1), axes = FALSE, xlab = "", ylab = "", col = myCol)
  Labels <- c("Min", "Max")
  axis(1, at = c(0,40), labels = Labels, las = 1)

  if(writeToFile == TRUE) {
    dev.off()
  }
  else {
    .wait()
  }
}

.plotCardStats = function(plotMat, plotTitle,  writeToFile, statType) 
{
  if (writeToFile == TRUE) {
    jpegTitle <- paste(plotTitle,".jpg",sep="")
    jpeg(jpegTitle)
  }
  myPal <- brewer.pal(5, "RdYlGn")
  myCol <- colorRampPalette(myPal[c(1:5,4:1)])(128)
  background <- "black"
  n <- nrow(plotMat)
  m <- ncol(plotMat)
  rname <- 1:n
  cname <- 1:m
  rname <- as.character(rname)
  cname <- as.character(cname)

  layout(matrix(c(1, 2), 2, 1), heights = c(4,1))
  op <- par(mar = c(0, 4, 4, 4))
  plot.new()
  plot.window(c(0, m), c(0, n), asp = 1)
  xlabwidth <- max(strwidth(rname, cex = 1))
  ylabwidth <- max(strwidth(cname, cex = 1))
  plot.window(c(-xlabwidth + 0.5, m + 0.5), c(0, n + 1 + ylabwidth), asp = 1, xlab="", ylab="")
  rect(0.5, 0.5, m + 0.5, n + 0.5, col = background)      #background color
  title(plotTitle)
  segments(rep(0.5, n + 1), 0.5 + 0:n, rep(m + 0.5, n + 1), 0.5 + 0:n, col = "gray")
  segments(0.5 + 0:m, rep(0.5, m + 1), 0.5 + 0:m, rep(n + 0.5, m), col = "gray")
  cols <- .computeColors(plotMat)
  bg <- myCol[cols]
  symbols(rep(1:m, each = n), rep(n:1, m), add = TRUE, inches = F, circles = rep(0.4, (m*n)), bg = as.vector(bg))
  
  x.bar <- seq(from = -3, to = 3, length = length(myCol))
  par(mar = c(5.1, 4.1, 1, 2))
  image(x.bar, 1, matrix(x.bar, length(x.bar), 1), axes = FALSE, xlab = "", ylab = "", col = myCol)
  x.small <- c(-3, -1.5, 0, 1.5, 3)
  if(statType == "parametric") Labels <- c("<=-3*SD", "1.5*SD", "mean", "1.5*SD", ">=3*SD")
  else  Labels <- c("<=-3*MAD", "1.5*MAD", "median", "1.5*MAD", ">=3*MAD")
  axis(1, at = x.small, labels = Labels, las = 1)

  if (writeToFile == TRUE) {
    dev.off()
  }
  else {
    .wait()
  }
}

.computeColors <- function(z){
  bound <- 3
  dispersion <- 1
  z0 <- z
  z[z0 < -bound] <- rep(1, sum(z0 < -bound))
  z[z0 >= -bound & z0 < -dispersion] <- abs(ceiling(abs(z0[z0 >= -bound & z0 < -dispersion] + dispersion)/(bound-dispersion)*32)-33)
  z[z0 >= -dispersion & z0 < 0] <- abs(ceiling(abs(z0[z0 >= -dispersion & z0 < 0])/dispersion*32)-33)+32
  z[z0 >= 0 & z0 < dispersion] <- ceiling(z0[z0 >= 0 & z0 < dispersion]/dispersion*32)+64
  z[z0 >= dispersion & z0 < bound] <- ceiling((z0[z0 >= dispersion & z0 <= bound] - dispersion)/(bound-dispersion)*32)+96
  z[z0 >= bound] <- rep(1, sum(z0 >= bound))
  z
}
.wait <- function() {
  par(ask=TRUE)
}
