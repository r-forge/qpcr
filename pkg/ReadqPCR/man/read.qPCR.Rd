\name{read.qPCR}
\alias{read.qPCR}
\alias{read_qPCRBatch}
\alias{checkValidqPCRFilename}
\title{ Read Taqman qPCR data into ExpressionSet }
\description{
  Reads RT-qPCR data in format specified in the ReadqPCR vignette to object of class
  \code{"ExpressionSet"}.
}
\usage{
read.qPCR(filename = character(0),
          phenoData = new("AnnotatedDataFrame"),
          notes = "",
          verbose = FALSE)
}
\arguments{
  \item{filename}{file name in a character vector.}
  \item{phenoData}{an \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
     object, a \code{character} of length one, or a \code{data.frame}.}
  \item{notes}{notes.}
  \item{verbose}{verbosity flag.}
}
\details{
  Reads a qPCR data file formatted to specifications in the ReadqPCR vignette
}
\value{Object of class \code{"ExpressionSet"}.}
%\references{ ~put references to the literature/web site here ~ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk} }
%\note{}
\seealso{ \code{\link[Biobase]{ExpressionSet-class}} }
%\examples{
%}
\keyword{classes}
