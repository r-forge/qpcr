\name{combineTechReps}
\alias{combineTechReps}
\alias{combineTechReps,qPCRBatch-method}
\title{ Takes expression set of qPCR values containing technical replicates and combines them. }
\description{
Takes expression set of qPCR values containing technical replicates and combines them.
}
\usage{
combineTechReps(qPCRBatch)
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data, read in by a ReadqPCR function and containing technical reps, denoted by _TechRep.n suffix.
}
}
\details{
  Takes expression set of qPCR values containing technical replicates and combines them.
}
\value{
 qPCRBatch with same number of samples, but with features reduced by the number of technical replicates that were present on the card.
}
%\references{ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
%\seealso{}
%\examples{
%}
\keyword{data}
