\name{makeAllNAs}
\alias{makeAllNAs}
\alias{makeAllNAs,qPCRBatch-method}
\title{ Make all NAs when number of NAs above a given threshold }
\description{ Make all NAs when number of NAs above a given threshold
}
\usage{
makeAllNAs(qPCRBatch, contrastM, sampleMaxM)
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data. 
}
  \item{contrastM}{ Contrast Matrix like that used in Limma. Columns represent the different samples types, rows are the different samples,
with a 1 or 0 in the matrix indicating which sample types the different samples belong to.
}
  \item{sampleMaxM}{ Sample Max Matrix. Columns represent the different sample types. There is one value per column, which represents the max number of NAs allowed for that sample type.
}
}
\details{
  Make all NAs when number of NAs above a given threshold
}
\value{
  qPCRBatch object with a new exprs slot, everything else equal
}
%\references{ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
%\seealso{}
%\examples{
%}
\keyword{data}
