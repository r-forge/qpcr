setGeneric("PseudoPlot",
  function(qPCRSet, plotType="Cts.Values", writeToFile=FALSE, cutOff=40)
    standardGeneric("PseudoPlot")
)
setMethod("PseudoPlot", signature = "qPCRSet", definition =
  function(qPCRSet, plotType, writeToFile, cutOff) {
    ctsMat <- exprs(qPCRSet)
    ctsMat[is.na(ctsMat)] <- cutOff # Cutoff point..Could change for different platforms
    orderMat <- exprs.well.order(qPCRSet)
    plateVec <- as.vector(gsub("-.*", "", orderMat))
    wellVec <- as.numeric(gsub(".*-", "", orderMat))
    minVal <- 0

    if (plotType == "Cts.Values") {
      maxVal <- round(max(ctsMat, na.rm=TRUE), 2)
      for (plate in unique(plateVec)) { # for each plate
        plotTitle <- paste(plotType, "for plate:", plate)
        orderedVals <- ctsMat[plateVec == plate][order(wellVec[plateVec == plate])]

        plotMat <- matrix(orderedVals, nrow=16, byrow=TRUE)
        .plotCard(plotMat, plotTitle, minVal, maxVal)
      }
    }
    else if (plotType == "Plate.Residuals") {
      for (plate in unique(plateVec)) {
        plotTitle <- paste(plotType, "for plate:", plate)
        orderedVals <- ctsMat[plateVec == plate][order(wellVec[plateVec == plate])]
        dispersion <- sd(orderedVals)
        dispersion <- sd(as.vector(orderedVals))
        meanSubbedVals <-  orderedVals - mean(orderedVals)
        plotMat <- matrix(meanSubbedVals, nrow=16, byrow=TRUE)
        .plotCard(plotMat, plotTitle, dispersion, writeToFile)
      }
    }
    else if (plotType == "Detector.Residuals") {
      valMat <- abs(ctsMat - apply(ctsMat, 1, mean, na.rm=TRUE)) # take the avg values from the Cts vals
      maxVal <- round(max(valMat, na.rm=TRUE), 2)
      for (plate in unique(plateVec)) {
        title <- paste(plotType, "for plate:", plate)
        orderedVals <- valMat[plateVec == plate][order(wellVec[plateVec == plate])]
        plotMat <- matrix(orderedVals, nrow=16, byrow=TRUE)
        .plotCard(plotMat, title, minVal, maxVal)
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
        maxVal <- round(max(plotMat, na.rm=TRUE), 2)
        .plotCard(plotMat, title, minVal, writeToFile)
      }
    }
  }
)

.plotCard = function(plotMat, plotTitle, dispersion, writeToFile) 
{
  if (writeToFile == TRUE) {
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
layout(matrix(c(1, 2), 2, 1), heights = c(4,1))
op <- par(mar = c(0, 4, 4, 4))
  plot.new()
cat("hereo")
  plot.window(c(0, m), c(0, n), asp = 1)
  xlabwidth <- max(strwidth(rname, cex = 1))
  ylabwidth <- max(strwidth(cname, cex = 1))
  plot.window(c(-xlabwidth + 0.5, m + 0.5), c(0, n + 1 + ylabwidth), asp = 1, xlab="", ylab="")
  rect(0.5, 0.5, m + 0.5, n + 0.5, col = background)      ##background color
  title(plotTitle)
  segments(rep(0.5, n + 1), 0.5 + 0:n, rep(m + 0.5, n + 1), 0.5 + 0:n, col = "gray")
  segments(0.5 + 0:m, rep(0.5, m + 1), 0.5 + 0:m, rep(n + 0.5, m), col = "gray")
  cols <- .computeColors(plotMat, dispersion)

  bg <- myCol1[cols]
print(cols)
print(bg)
  symbols(rep(1:m, each = n), rep(n:1, m), add = TRUE, inches = F, circles = rep(0.4, (m*n)), bg = as.vector(bg))
x.bar <- seq(from = -3*dispersion, to = 3*dispersion, length = length(myCol1))
par(mar = c(5.1, 4.1, 1, 2))
image(x.bar, 1, matrix(x.bar, length(x.bar), 1), axes = FALSE, xlab = "", ylab = "", col = myCol1)
x.small <- c(-3*dispersion, -1.5*dispersion, 0, 1.5*dispersion, 3*dispersion)
Labels <- c("<=-3*SD", "SD", "mean", "SD", ">=3*SD")
axis(1, at = x.small, labels = Labels, las = 1)
cat("cols:\n")

  if (writeToFile == TRUE) {
    dev.off()
  }
  else {
    .wait()
  }
}

library(RColorBrewer)
myPal <- brewer.pal(5, "RdYlGn")
myCol1 <- colorRampPalette(myPal[c(1:5,4:1)])(128)

.computeColors <- function(z, dispersion){
    bound <- dispersion*3
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
