\name{normByHKG}
\alias{normByHKG}
\alias{normaliseByHKG}
\title{ Normalise qPCR data with housekeeping genes }
\description{
  Normalise qPCR eset using housekeeping genes as control.
}
\usage{
normByHKG(qPCRSet, hkg, verbose = FALSE)
}
\arguments{
  \item{qPCRSet}{ Expression set containing qPCR data. 
}
  \item{hkg}{ String containing the name of th name of the housekeeping gene which will be used to normalise the rest of the genes.
}
  \item{verbose}{ Verbosity flag.}
}
\details{
  Takes expression set of qPCR values and normalises them using different housekeeping genes. Returns seperate sets of values for each housekeeping gene given.
}
\value{
  Object of class \code{"qPCRSet"}
}
%\references{ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
%\seealso{}
%\examples{
%}
\keyword{data}
