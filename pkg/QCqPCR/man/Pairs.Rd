\name{qPCRPairs}
\alias{qPCRPairs}
\alias{qPCRPairs,qPCRBatch-method}
\title{ Make pairwise graphs for qPCR data from qPCR eset }
\description{
  Takes set of qPCR data and makes pairwise graphs between all the samples.
}
\usage{
qPCRPairs(qPCRBatch, writeToFile = FALSE)
}
\arguments{
    \item{qPCRBatch}{ Expression set containing qPCR data.}
    \item{writeToFile}{ Option to write to file}
}
\details{
  reads in qPCR data and pairwise graphs, optionally saving them to file
}
\section{Methods}{
  \describe{
    \item{qPCRPairs}{(qPCRBatch): ... }
  }
}
\value{produces graph and optionally writes to ping if called.}
%\references{ ~put references to the literature/web site here ~ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk} }
%\examples{
%}
\keyword{classes}
