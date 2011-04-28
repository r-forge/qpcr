\name{plotVsHkg}
\alias{plotVsHkg}
\alias{plotVsHkg,qPCRBatch-method}
\title{ Plot Ct Values when normalised against different housekeeping genes }
\description{
Takes normalised qPCR data and plot the different Ct values from the different house keeping genes
}
\usage{
plotVsHkg(qPCRBatch, \dots)

\S4method{plotVsHkg}{qPCRBatch}(qPCRBatch, hkgs, writeToFile = FALSE) 
}
\arguments{
  \item{qPCRBatch}{ data set containing data on different Ct values normed to the different HKGs }
  \item{\dots}{ Extra arguments, detailed below }
  \item{hkgs}{ Vector of housekeeping genes. Must all be present in the exprs matrix } 
  \item{writeToFile}{ Option to write output to file }
}
\details{
  Takes normalised qPCR data and plot the different Ct values from the different house keeping genes.
  Plots HKGs against each other on the same graph, ordered by a single HKG
}
\value{produces graph and optionally writes to png if called.}
%\references{ ~put references to the literature/web site here ~ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk} }
%\examples{
%}
\keyword{classes}
