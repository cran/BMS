\name{fullmodel.ssq}
\alias{fullmodel.ssq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ OLS Statistics for the Full Model Including All Potential Covariates }
\description{
  A utility function for reference: Returns a list with R2 and sum of squares for the OLS model encompassing all potential covariates that are included in a bma object.
}
\usage{
fullmodel.ssq(yX.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yX.data}{ a bma object (cf. \code{\link{bms}}) - alternatively a \link{data.frame} or \link{matrix} whose first column is the dependent variable} 
}
\value{
  Returns a list with some basic OLS statistics
  \item{R2}{The R-squared of the full model}
  \item{ymy}{The sum of squares of residuals of the full model}
  \item{ypy}{The explained sum of squares of the full model}
  \item{yty}{The sum of squares of the (demeaned) dependent variable}
  \item{Fstat}{The F-statistic of the full model}
}
\author{ Martin Feldkircher and Stefan Zeugner }
\note{ This function is just for quick comparison; for proper OLS estimation consider \code{\link{lm}} }
\seealso{ \code{\link{bms}} for creating bma objects, \code{\link{lm}} for OLS estimation 

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
data(datafls)
mm=bms(datafls)

fullmodel.ssq(mm)

#equivalent:
fullmodel.ssq(mm$X.data)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
