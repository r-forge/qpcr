\name{pairsqPCR}
\alias{pairsqPCR}
\title{ Make pairwise graphs for qPCR data from qPCR eset }
\description{
  Takes set of qPCR data and makes pairwise graphs between all the samples.
}
\usage{
pairsqPCR(qPCRSet, saveToFile = FALSE)
}
\arguments{
    \item{qPCRSet}{ Expression set containing qPCR data.}
    \item{saveToFile} { Option to save to file}
}
\details{
  reads in qPCR data and pairwise graphs, optionally saving them to file
}
\value{produces graph and optionally writes to ping if called.}
%\references{ ~put references to the literature/web site here ~ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk} }
%\examples{
%}
\keyword{classes}
