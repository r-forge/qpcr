name{read.taqmanbatch}
\alias{read.taqmanbatch}
\alias{read_TaqBatch}
\alias{checkValidTaqmanFilenames}
\title{ Read Taqman qPCR data into ExpressionSet }
\description{
  Reads Taqman RT-qPCR data to object of class
  \code{"ExpressionSet"}.
}
\usage{
read.taqmanbatch(\dots,
                 filenames = character(0),
                 phenoData = new("AnnotatedDataFrame"),
                 notes = "", 
                 verbose = FALSE)
}
\arguments{
  \item{\dots}{file names separated by comma.}
  \item{filenames}{file names in a character vector.}
  \item{phenoData}{an \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
     object, a \code{character} of length one, or a \code{data.frame}.}
  \item{notes}{notes.}
  \item{verbose}{verbosity flag.}
}
\details{
  Reads a single Taqman data file.
}
\value{Object of class \code{"ExpressionSet"}.}
%\references{ ~put references to the literature/web site here ~ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk} }
%\note{}
\seealso{ \code{\link[Biobase]{ExpressionSet-class}} }
%\examples{
%}
\keyword{classes}

