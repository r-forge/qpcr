\name{deltaDeltaCt}
\alias{deltaDeltaCt}
\alias{deltaDeltaCt-methods}
\alias{deltaDeltaCt,qPCRBatch-method}
\title{ Perform normalization and differential expression with given housekeeping gene }
\description{
  Normalise qPCR eset using housekeeping genes as control, then perform differential expression analysis 
  using the   delta delta Ct method.
  Differs from the deltaDeltaAvgCt method in that the hkg values are subtracted and the sd is 
  calculated on the differences (i.e. paired).
  Suitable when housekeeping genes are from same wells/sample as the other detectors
}
\usage{
  deltaDeltaCt(qPCRBatch,\dots)

  \S4method{deltaDeltaCt}{qPCRBatch}(qPCRBatch, maxNACase=0, maxNAControl=0, hkgs, contrastM, 
             case, control, paired=TRUE, hkgCalc="arith", statCalc="arith")
}
\arguments{
  \item{qPCRBatch}{ qPCR-specific expression set, containing qPCR data. }
  \item{\dots}{ Extra arguments, detailed below }
  \item{maxNACase}{ Maximum number of NA values allowed before a detector's reading is discarded for samples designated as case. }
  \item{maxNAControl}{ Maximum number of NA values allowed before a detector's reading is discarded for samples designated as control. }
  \item{hkgs}{ String containing the name of th name of the housekeeping gene which will be used to normalise the rest of the genes. }
  \item{contrastM}{ A binary matrix which designates case and control samples. }
  \item{case}{ The name of the column in contrastM that corresponds to the case samples. }
  \item{control}{ The name of the column in contrastM that corresponds to the control samples. }
  \item{paired}{ Logical - if TRUE the detectors and housekeepers in the same sample will be paired for calculating standard deviation, effectively meaning we will be calculating standard deviation of the differences. If FALSE, there will be no pairing, and standard deviation will be pooled between the detector and housekeepers. }
  \item{hkgCalc}{ String - either "arith" or "geom", details how the different housekeeper genes should be combined - either by using the arithmetic or geometric mean. }
  \item{statCalc}{ String - either "arith" or "geom", details how genes should be combined - either by using the arithmetic or geometric mean. }
}
\details{
  Takes expression set of qPCR values and normalises them using different housekeeping genes. Returns seperate sets of values for each housekeeping gene given.
}
\value{
  matrix with columns containing the detector ids, 2^delta Ct values for the sample of interest and the callibrator sample, alongside their respective standard deviations, the 2^delta delta Ct values and the minimum and maximum values (ie the values that are 1 sd away )
}
\references{
  Kenneth Livak, Thomase Schmittgen (2001).
  Analysis of Relative Gene Expression Data Using Real-Time Quantitative
  PCR and the 2^DDCt Method.
  Methods 25, 402-408, 2001
  \url{http://www.ncbi.nlm.nih.gov/pubmed/11846609}
}
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
\seealso{\code{selectHKs}, \code{deltaCt}}
\examples{
  path <- system.file("exData", package = "NormqPCR")
  taqman.example <- file.path(path, "example.txt")
  qPCRBatch.taqman <- read.taqman(taqman.example)
  hkg <- "Actb-Rn00667869_m1"

  contM <- cbind(c(0,0,1,1,0,0,1,1),c(1,1,0,0,1,1,0,0))
  colnames(contM) <- c("interestingPhenotype","wildTypePhenotype")
  rownames(contM) <- sampleNames(qPCRBatch.taqman)

  ddCt.taqman <- deltaDeltaCt(qPCRBatch = qPCRBatch.taqman, maxNACase=1, maxNAControl=1, hkg=hkg, contrastM=contM, case="interestingPhenotype", control="wildTypePhenotype", statCalc="geom", hkgCalc="arith")
  head(ddCt.taqman)
}
\keyword{data}
