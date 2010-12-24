\name{combineTechReps}
\alias{combineTechReps}
\alias{combineTechReps,qPCRBatch-method}
\title{ Combines Technical Replicates }
\description{
Takes expression set of qPCR values containing technical replicates and combines them.
}
\usage{
combineTechReps(qPCRBatch)
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data, read in by a ReadqPCR function and containing technical reps, denoted by \code{_TechRep.n} suffix.
}
}
\details{
  Takes \code{exprs} of qPCR values containing technical replicates and combines them.
}
\value{
 \code{qPCRBatch} with same number of samples, but with less features, since all technical replicates are replaced with a single value of their means.
}
%\references{ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
%\seealso{}
\examples{
  path <- system.file("exData", package = "NormqPCR")
  qPCR.example.techReps <- paste(path,"/qPCR.techReps.txt", sep = "")
  qPCRBatch.qPCR.techReps <- read.qPCR(qPCR.example.techReps)
  rownames(exprs(qPCRBatch.qPCR.techReps))
  combinedTechReps <- combineTechReps(qPCRBatch.qPCR.techReps)
  rownames(exprs(combinedTechReps))
}
\keyword{data}
