\name{is.bma}
\alias{is.bma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Tests for a 'bma' Object }
\description{
  tests for objects of class "bma"
}
\usage{
is.bma(bmao)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bmao}{ a 'bma' object: see 'value'  }
}
%\details{~~ If necessary, more details than the description above ~~}
\value{
  Returns \code{TRUE} if bmao is of class 'bma', \code{FALSE} otherwise.  
}

\author{ Martin Feldkircher and Stefan Zeugner}
\seealso{ 'Output' in \code{\link{bms}} for the structure of a 'bma' object 

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
 data(datafls)
 mm=bms(datafls,burn=1000, iter=4000)
 is.bma(mm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
