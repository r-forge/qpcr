\name{taq2eset}
\alias{taq2eset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
taq2eset(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ ~~Describe \code{filename} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(filename) {
raw = read.delim(filename)
levels(raw$Sample) = gsub(" ","_",levels(raw$Sample))
samples = levels(raw$Sample)
detectors = levels(raw$Detector)

eset = data.frame(detectors,row.names=1)
for (sample in samples) {
	Cts = data.frame(raw$Detector[raw$Sample == sample],as.numeric(as.character(raw$Ct[raw$Sample == sample])),row.names=1)
	eset = data.frame(merge(eset,Cts,by="row.names"),row.names=1)
  }
names(eset) = samples
eset = new("ExpressionSet",exprs=as.matrix(eset)) # make it an expression set

return(eset)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
