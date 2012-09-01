\name{variable.names.zlm}
\alias{variable.names.zlm}
\alias{vcov.zlm}
\alias{logLik.zlm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Variable names and design matrix }
\description{
  Simple utilities retrieving variable names and design matrix from a bma object
}
\usage{

 \method{variable.names}{zlm}(object,...)
 \method{vcov}{zlm}(object, include.const = FALSE, ...)
 \method{logLik}{zlm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{bma} object (as produced by \code{\link{bms}}) }
  \item{include.const}{ Whether the variance-covariance matrix returned by \code{vcov.zlm} should also include a line and row for the intercept (which will be NA for most priors) }
  \item{\ldots}{further arguments passed to or from other methods }
}
\details{
 \code{variable.names.zlm}: method \code{\link{variable.names}} for a \code{\link{zlm}} model. \cr
 \code{vcov.zlm}: the posterior variance-covariance matrix of the coefficients of a \code{\link{zlm}} model - cf. \code{\link{vcov}} \cr
 \code{logLik.zlm}: a \code{\link{zlm}} model's log-likelihood \code{p(y|M)} according to the implementation of the respective coefficent prior \cr
}

\author{ Martin Feldkircher and Stefan Zeugner }
\seealso{ \code{\link{zlm}} for creating \code{zlm} objects 

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
 data(datafls)
  
 zz=zlm(datafls)
 variable.names(zz)
 vcov(zz)
 logLik(zz)
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
