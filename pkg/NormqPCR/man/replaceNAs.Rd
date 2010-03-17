\name{replaceNAs}
\alias{replaceNAs}
\alias{replaceNAs,qPCRBatch-method}
\title{ Replace NAs with a given value }
\description{ Replace NAs with a given value
}
\usage{
replaceNAs(qPCRBatch, newNA)
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data. 
}
  \item{newNA}{ The new value to replace the NAs with
}
}
\details{
  Replaces NA values in the exprs slot of the qPCRBatch object with a given number
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
