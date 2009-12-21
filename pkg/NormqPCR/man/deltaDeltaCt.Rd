\name{deltaDeltaCt}
\alias{deltaDeltaCt}
\alias{deltaDeltaCt,qPCRBatch-method}
\title{ Perform normalization and differential expression with given housekeeping gene }
\description{
Normalise qPCR eset using housekeeping genes as control, then perform differential expression analysis using the delta delta Ct method.
Differs from the deltaDeltaAvgCt method in that the hkg values are subtracted and the sd is calculated on the differences (i.e. paired).
Suitable when housekeeping genes are from same wells/sample as the other detectors
}
\usage{
deltaDeltaCt(qPCRBatch, maxNACase=0, maxNAControl=0, hkg, contrastM, case, control)
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data. 
}
  \item{maxNACase}{ Maximum number of NA values allowed before a detector's reading is discarded for samples designated as case
}
  \item{maxNAControl}{ Maximum number of NA values allowed before a detector's reading is discarded for samples designated as control
}
  \item{hkg}{ String containing the name of th name of the housekeeping gene which will be used to normalise the rest of the genes.
}
  \item{contrastM}{ A binary matrix which designates case and control samples
}
  \item{case}{ The name of the column in contrastM that corresponds to the case samples
}
  \item{control}{ The name of the column in contrastM that corresponds to the control samples
}
}
\details{
  Takes expression set of qPCR values and normalises them using different housekeeping genes. Returns seperate sets of values for each housekeeping gene given.
}
\value{
  matrix with columns containing the detector ids, the 2^delta delta Ct values and the minimum and maximum values (ie the values that are 1 sd away )
}
%\references{ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
%\seealso{}
%\examples{
%}
\keyword{data}
