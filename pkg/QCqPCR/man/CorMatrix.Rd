\name{qPCRCorMatrix}
\alias{qPCRCorMatrix}
\alias{qPCRCorMatrix,qPCRBatch-method}
\title{ Make pairwise matrix of correlation values between samples in qPCR eset }
\description{
Takes set of qPCR data and makes matrix of correlation values using pearson's coefficient
}
\usage{
  qPCRCorMatrix(qPCRBatch)
  \S4method{qPCRCorMatrix}{qPCRBatch}(qPCRBatch) 
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data.}
}
\details{
  reads in qPCR data and returns  matrix of correlation values using Pearon's
  coefficient, 
}
\section{Methods}{
  \describe{
    \item{qPCRCorMatrix}{(qPCRBatch): ... }
  }
}
\value{produces nxn matrix of correlation values.}
%\references{ ~put references to the literature/web site here ~ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk} }
%\examples{
%}
\keyword{classes}
