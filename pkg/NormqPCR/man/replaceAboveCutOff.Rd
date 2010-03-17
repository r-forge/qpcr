\name{replaceAboveCutOff}
\alias{replaceAboveCutOff}
\alias{replaceAboveCutOff,qPCRBatch-method}
\title{ Replace values above a given threshold with a new value }
\description{ Replace values above a given threshold with a new value
}
\usage{
replaceAboveCutOff(qPCRBatch, newVal=NA, cutOff=38)
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data. 
}
  \item{newVal}{ The new value with which to replace the values above the cutoff
}
  \item{cutOff}{ the minimal threshold above which the values will be replaced
}
}
\details{
  Replaces values in the exprs slot of the qPCRBatch object that are above a threshold value with a new number
}
\value{
  qPCRBatch object with a new exprs slot
}
%\references{ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
%\seealso{}
%\examples{
%}
\keyword{data}
