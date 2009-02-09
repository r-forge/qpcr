\name{pairsqPCR}
\alias{pairsqPCR}
\title{ Make pairwise graphs for qPCR data from qPCR eset }
\description{
  Takes eset of qPCR data and makes pairwise graphs between all the samples.
}
\usage{
pairsqPCR(eset, saveToFile = FALSE)
}
\arguments{
    \item{eset}{ Expression set containing qPCR data.}
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
