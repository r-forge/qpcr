\name{NormqPCR-package}
\alias{NormqPCR-package}
\alias{NormqPCR}
\docType{package}
\title{
  Functions for normalisation of RT-qPCR data.
}
\description{
Functions for normalisation of real-time quantitative PCR data.
}
\details{
\tabular{ll}{
Package: \tab NormqPCR\cr
Type: \tab Package\cr
Version: \tab 0.99.3\cr
Date: \tab 2011-07-15\cr
Depends: \tab R(>= 2.10.0), stats, RColorBrewer, Biobase, methods, ReadqPCR \cr
License: \tab LGPL-3\cr
LazyLoad: \tab yes \cr
LazyData: \tab yes \cr

}
require(NormqPCR)
}
\author{ Matthias Kohl, James Perkins

Maintainer: Matthias Kohl  \email{Matthias.Kohl@stamats.de}}
\references{
  Jo Vandesompele, Katleen De Preter, Filip Pattyn et al. (2002). Accurate 
  normalization of real-time quantitative RT-PCR data by geometric averaging 
  of multiple internal control genes. 
  Genome Biology 2002. 3(7):research0034.1-0034.11.
  \url{http://genomebiology.com/2002/3/7/research/0034/}
  
  Claus Lindbjerg Andersen, Jens Ledet Jensen and Torben Falck Orntoft (2004).
  Normalization of Real-Time Quantitative Reverse Transcription-PCR Data: 
  A Model-Based Variance Estimation Approach to Identify Genes Suited for
  Normalization, Applied to Bladder and Colon Cancer Data Sets. 
  CANCER RESEARCH 64, 5245-5250, August 1, 2004.
  \url{http://cancerres.aacrjournals.org/cgi/content/full/64/15/5245}

  Kenneth Livak, Thomase Schmittgen (2001).
  Analysis of Relative Gene Expression 
  Data Using Real-Time Quantitative PCR and the 2^ddCt Method.
  Methods 25, 402-408, 2001
  \url{http://www.ncbi.nlm.nih.gov/pubmed/11846609}
}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
## some examples are given in the vignette
\dontrun{
library(NormqPCR)
vignette("NormqPCR")
}
}
\keyword{package}
