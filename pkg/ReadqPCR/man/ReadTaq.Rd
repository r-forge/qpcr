\name{ReadTaq}
\alias{ReadTaq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ReadTaq(filename = character())
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
function(filename = character()) {

# If no filename given look for .taq files in current directory
if (length(filename) == 0) {
	taqfile.path <- getwd()
	files = list.files(taqfile.path)
	filename = files[grep("\\.[tT][aA][qQ]$",files)]
	if (length(filename)>1) {
		stop("Too many Taq files in current directory","\n")
	}
  }
if (length(filename) == 0) {
	stop("No Taqman filenname specified and no Taqman files in current directory:","\n")
  }
if (length(filename)>1) {
	stop ("Too may Taqman files specified")
  }
taq2eset(filename)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
