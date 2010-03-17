\name{deltaCt}
\alias{deltaCt}
\alias{deltaCt,qPCRBatch-method}
\title{ Perform normalization with a given housekeeping gene }
\description{
Normalise qPCR eset using a given housekeeping gene as control, then perform differential expression analysis using the delta delta Ct method
}
\usage{
deltaCt(qPCRBatch, hkgs)
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data. 
}
  \item{hkgs}{ String containing the name of the name of the housekeeping gene which will be used to normalise the rest of the genes.
}
}
\details{
  Takes expression set of qPCR values and normalises them using a housekeeping gene. Returns a qPCRBatch with exprs set of the same dimensions but with the given hkg value subtracted.
}
\value{
 qPCRBatch with exprs set of the same dimensions but with the given hkg value subtracted.

}
%\references{ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
%\seealso{}
%\examples{
%}
\keyword{data}
