\name{produceHkgsDF}
\alias{produceHkgsDF}
\title{ Produce Matrix Comparing the Relative Performance of the Different Houskeeping Genes. }
\description{
  Produce Matrix Comparing the Relative Performance of the Different Houskeeping Genes
  This servers as the input for the plotVsHkg QC function.
}
\usage{
produceHkgsDF(qPCRSet, hkgs, design, cutoff = 40, verbose = FALSE)
}
\arguments{
  \item{qPCRSet}{ Expression set containing qPCR data. 
  }
  \item{hkgs}{ String containing the name of th name of the housekeeping gene which will be used to normalise the rest of the genes.
  }
  \item{design}{ Vector specifying which of the columns of the exprs represent case and which represent control.
  }
  \item{cutoff}{ gives a cutoff of cycles after which it can be assumed that there is nothing present.
  }
  \item{verbose}{ Verbosity flag.
  }
}
\details{
  Produce Matrix Comparing the Relative Performance of the Different Houskeeping Genes.
  This servers as the input for the plotVsHkg QC function.
}
\value{
  Object of class \code{"qPCRSet"}
}
%\references{ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
%\seealso{}
%\examples{
%}
\keyword{data}
