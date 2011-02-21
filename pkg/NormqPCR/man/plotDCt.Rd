\name{plotDCt}
\alias{plotDCt}
\alias{plotDCt,qPCRBatch-method}
\title{Plot delta Ct values}
\description{
Takes expression set of qPCR values and plots Delta Ct after subtraction of normalisation factor (NF) or housekeeper gene
}
\usage{
plotDCt(..., ddCtTable, detectors="", statCalc="arith")
}
\arguments{
  \item{\dots}{ Additional arguments to pass to the plot command.
}
  \item{ddCtTable}{ output of deltaDeltaCt.
}
  \item{detectors}{ The detectors to show, if left blank will show all detectors.
}
  \item{statCalc}{ The type of statistic to show.
}

}
\details{
  Takes ddCtTable output from the deltaDeltaCt() command and plots the Ct values plus variance for a given number of detectors.
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
