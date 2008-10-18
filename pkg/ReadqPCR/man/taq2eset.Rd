\name{taq2eset}
\alias{taq2eset}
\title{ Convert Taqman data to ExpressionSet }
\description{
  Function to convert raw Taqman data to object of class 
  \code{"ExpressionSet"}.
}
\usage{
taq2eset(filename)
}
\arguments{
  \item{filename}{ character, name of Taqman data file.}
}
\details{
  Convert a single Taqman RT-qPCR file to object of class
  \code{"ExpressionSet"}.
}
\value{Object of class \code{"ExpressionSet"}.}
%\references{ ~put references to the literature/web site here ~ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk} }
%\note{}
\seealso{ \code{\link[Biobase]{ExpressionSet-class}} }
%\examples{
%}
\keyword{classes}
