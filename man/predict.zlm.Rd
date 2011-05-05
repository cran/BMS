\name{predict.zlm}
\alias{predict.zlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Predict Method for zlm Linear Model}
\description{
  Expected value (And standard errors) of predictions based on 'zlm' linear Bayesian model under Zellner's g prior
}
\usage{
## S3 method for class 'zlm':
 \method{predict}{zlm}(object, newdata = NULL, se.fit = FALSE, ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ a zlm linear model object - see \code{\link{zlm}} }
  \item{newdata}{An optional data.frame, matrix or vector containing variables with which to predict. If omitted, then (the expected values of) the fitted values are returned.}
  \item{se.fit}{A switch indicating if the standard deviations for the predicted varaibles are required.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  A vector with (expected values of) fitted values.\cr
  If \code{se.fit} is \code{TRUE}, then the output is a list with the following elements:
  \item{fit}{ a vector with the expected values of fitted values}
  \item{std.err}{ a vector with the standard deviations of fitted values}
  \item{se.fit}{ a vector with the standard errors without the residual scale akin to \code{se.fit} in \code{\link{predict.lm}} }
  \item{residual.scale}{ The part from the standard deviations that involves the identity matrix. Note that \code{sqrt(se.fit^2+residual.scale^2)} yields \code{std.err}. }
  
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Martin Feldkircher and Stefan Zeugner }
\seealso{ \code{\link{bms}} for creating zlm objects, \code{\link{predict.lm}} for a comparable function, \code{\link{predict.bma}} for predicting with bma objects

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
 data(datafls)
 mm=zlm(datafls,g="EBL")
 
 predict(mm) #fitted values 
 predict(mm, newdata=1:41) #prediction based on a 'new data point' 
 
 #prediction based on a 'new data point', with 'standard errors'
 predict(mm, newdata=datafls[1,], se.fit=TRUE) 
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
