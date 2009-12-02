\name{plotVsHkg}
\alias{plotVsHkg}
\title{ Plot Ct Values when normalised against different housekeeping genes }
\description{
  Takes normalised qPCR data and plot the different Ct values from the different house keeping genes
}
\usage{
plotVsHkg(qPCRSet, saveToFile = FALSE)
}
\arguments{
    \item{qPCRSet}{ data set containing data on different Ct values normed to the different HKGs }
    \item{saveToFile}{ Option to save to file }
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
