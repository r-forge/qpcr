\name{plotDdCt}
\alias{plotDdCt}
\alias{plotDdCt,qPCRBatch-method}
\title{Plot delta Ct values}
\description{
Takes expression set of qPCR values and plots Delta delta Ct after subtraction of normalisation factor (NF) or housekeeper gene and normalisation by subtraction of control
}
\usage{
plotDdCt(..., ddCtTable, detectors="", logFC=TRUE)
}
\arguments{
  \item{\dots}{ Additional arguments to pass to the plot command.
}
  \item{ddCtTable}{ output of deltaDeltaCt.
}
  \item{detectors}{ The detectors to show, if left blank will show all detectors.
}
  \item{logFC}{ plot fold change or log fold change.
}

}
\details{
  Takes ddCtTable output from the deltaDeltaCt() command and plots the ddCt values.
}
\value{
 Produces a barplot using barplot2
}
%\references{ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
%\seealso{}
%\examples{}
\keyword{plot,barplot,barplot2,barchart}
