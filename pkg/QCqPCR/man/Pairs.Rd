\name{qPCRPairs}
\alias{qPCRPairs}
\alias{qPCRPairs,qPCRBatch-method}
\title{ Make pairwise graphs for qPCR data from qPCR eset }
\description{
  Takes set of qPCR data and makes pairwise graphs between all the samples.
}
\usage{
qPCRPairs(qPCRBatch, plotType = "Sample", writeToFile = FALSE, pairsToPlot = "All")
}
\arguments{
    \item{qPCRBatch}{ Expression set containing qPCR data.}
    \item{plotType}{ type of plot. So far, \code{"Sample"} and \code{"Plate"}
      are implemented. Default is \code{"Sample"}. }
    \item{writeToFile}{ logical, write to jpeg-file. Default is \code{FALSE}. }
    \item{pairsToPlot}{ which pairs shall be plotted. Default is \code{"All"}. }
}
\details{
  generate pairwise graphs, optionally saving them to file.
}
\section{Methods}{
  \describe{
    \item{qPCRPairs}{(qPCRBatch): method for object of class \code{"qPCRBatch"}. }
  }
}
\value{produces graph and optionally writes to jpeg if called.}
%\references{ ~put references to the literature/web site here ~ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk} }
%\examples{
%}
\keyword{classes}
