\name{predict.zlm}
\alias{predict.zlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Predict Method for zlm Linear Model}
\description{
  Expected value of predictions based on 'zlm' linear Bayesian model under Zellner's g prior
}
\usage{
 \method{predict}{zlm}(object, newdata = NULL, ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ a zlm linear model object - see \code{\link{zlm}} }
  \item{newdata}{An optional data.frame, matrix or vector containing variables with which to predict. If omitted, then (the expected values of) the fitted values are returned.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  A vector with (expected values of) fitted values.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Martin Feldkircher and Stefan Zeugner }
\seealso{ \code{\link{bms}} for creating zlm objects, \code{\link{predict.lm}} for a comparable function

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
 data(datafls)
 mm=zlm(datafls,g="EBL")
 
 predict(mm) #fitted values 
 predict(mm, newdata=1:41) #prediction based on a 'new data point' 
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
