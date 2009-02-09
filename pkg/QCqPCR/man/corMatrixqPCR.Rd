\name{corMatrixqPCR}
\alias{corMatrixqPCR}
\title{ Make pairwise matrix of correlation values between samples in qPCR eset }
\description{
  Takes eset of qPCR data and makes matrix of correlation values using pearson's coefficient
}
\usage{
  corMatrixqPCR(eset)
}
\arguments{
    \item{eset}{ Expression set containing qPCR data.}
}
\details{
  reads in qPCR data and returns  matrix of correlation values using Pearon's
  coefficient, 
}
\value{produces nxn matrix of correlation values.}
%\references{ ~put references to the literature/web site here ~ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk} }
%\examples{
%}
\keyword{classes}
