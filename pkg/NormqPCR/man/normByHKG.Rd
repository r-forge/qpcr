\name{normByHKG}
\alias{normByHKG}
\alias{normaliseByHKG}
\title{ Normalise qPCR data with housekeeping genes }
\description{
  Normalise qPCR eset using housekeeping genes as control.
}
\usage{
normaliseByHKG(eset,hkgs,design,verbose = FALSE)
}
\arguments{
  \item{eset}{ Expression set containing qPCR data. 
}
  \item{hkgs}{ Character vector of housekeeping genes which will be used to normalise the rest of the genes, assuming there are at least 2 vectors.
}
  \item{design}{ Currently vector specifying which of the columns of the exprs represent case and which represent control. 
}
  \item{verbose}{ Verbosity flag.}
}
\details{
  Takes expression set of qPCR values and normalises them using different housekeeping genes. Returns seperate sets of values for each housekeeping  gene given.
}
\value{
  Object of class \code{"NormedqPCRSet"}
}
%\references{ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
%\seealso{}
%\examples{
%}
\keyword{data}
