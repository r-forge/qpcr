\name{PseudoPlot}
\alias{PseudoPlot}
\alias{PseudoPlot,qPCRBatch-method}
\title{ Plot Ct Values when normalised against different housekeeping genes }
\description{
Takes normalised qPCR data and plot the different Ct values from the different house keeping genes
}
\usage{
PseudoPlot(qPCRBatch, \dots)

\S4method{PseudoPlot}{qPCRBatch}(qPCRBatch, plotType = "Cts.Values", writeToFile = FALSE, 
           cutOff = NA, statType = "parametric", plateToPlot = "AllPlates") 
}
\arguments{
  \item{qPCRBatch}{ data set containing data on different Ct values normed to the different HKGs }
  \item{\dots}{ Extra arguments, detailed below }
  \item{plotType}{ yet to fill }
  \item{writeToFile}{ yet to fill }
  \item{cutOff}{ yet to fill }
  \item{statType}{ yet to fill }
  \item{plateToPlot}{ yet to fill }
}
\details{
  Takes normalised qPCR data and plot the different Ct values from the different house keeping genes.
  Plots HKGs against each other on the same graph, ordered by a single HKG
}
\section{Methods}{
  \describe{
    \item{PseudoPlot}{(qPCRBatch): method for object of class \code{"qPCRBatch"} }
  }
}

\value{produces graph and optionally writes to png if called.}
%\references{ ~put references to the literature/web site here ~ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk} }
%\examples{
%}
\keyword{classes}
