\name{print.bma}
\alias{print.bma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Printing bma Objects }
\description{
  Print method for objects of class 'bma'
}
\usage{
 \method{print}{bma}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a bma object - see \code{\link{bms}} }
  \item{\dots}{ additional arguments to \code{\link{print}} }
}
\value{
  prints calls to \code{coef.bma{x}} and \code{summary.bma(x)}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Martin Feldkircher and Stefan Zeugner }
\seealso{ \code{\link{coef.bma}} and \code{\link{summary.bma}} for the underlying functions, \code{\link{bms}} for creating bma objects, \code{\link{print}} 

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
 data(datafls)
 mm=bms(datafls,user.int=FALSE)
 print(mm)
 
 ## similar:
 coef(mm)
 summary(mm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
