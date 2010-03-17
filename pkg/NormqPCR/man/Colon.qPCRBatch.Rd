\name{Colon.qPCRBatch}
\alias{Colon.qPCRBatch}
\docType{data}
\title{ Data set of Andersen et al (2004) }
\description{
  This data set was used in Andersen et al (2004) to demonstrate normalization 
  of real-time quantitative RT-PCR data by geometric averaging of housekeeping genes.
}
\usage{data(geNorm)}
\format{
  A qPCRBatch object (see package ReadqPCR for full details) which extends expressionSet for RT-qPCR data
}
\details{
  The genes included in this data set were selected by screening 161 colon 
  sample expression profiles.
}
\source{
  The data set was obtained from
  \url{http://www.mdl.dk/Publications_sup1.htm}
  and then parsed into a different format and read into a qPCRBatch object using ReadqPCR packaged
}
\references{
  Claus Lindbjerg Andersen, Jens Ledet Jensen and Torben Falck Orntoft (2004).
  Normalization of Real-Time Quantitative Reverse Transcription-PCR Data: 
  A Model-Based Variance Estimation Approach to Identify Genes Suited for
  Normalization, Applied to Bladder and Colon Cancer Data Sets. 
  CANCER RESEARCH 64, 5245-5250, August 1, 2004.
  \url{http://cancerres.aacrjournals.org/cgi/content/full/64/15/5245}
}
\examples{
data(Colon.qPCRBatch)
str(Colon.qPCRBatch)
}
\keyword{datasets}
