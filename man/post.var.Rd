\name{post.var}
\alias{post.var}
\alias{post.pr2}
\alias{deviance.bma}
\alias{deviance.zlm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Posterior Variance and Deviance }
\description{
  Returns posterior residual variance, deviance, or pseudo R-squared, according to the chosen prior structure
}
\usage{
 post.var(object, exact = FALSE)
 post.pr2(object, exact = FALSE)
 \method{deviance}{bma}(object, exact = FALSE, ...)
 \method{deviance}{zlm}(object, ...)
 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{bma} object (as produced by \code{\link{bms}}) or a \code{\link{zlm}} object.}
  \item{exact}{ When \code{exact=FALSE}, then \code{deviance} will be based on MCMC frequencies, if \code{exact=TRUE} then it will be based on\cr 
    analytical posterior model probabilities - cf. argument \code{exact} in \code{\link{coef.bma}}.}
  \item{\ldots}{further arguments passed to or from other methods }
}
\details{
 \code{post.var}: Posterior residual variance as according to the prior definitions contained in \code{object} \cr
 \code{post.pr2}: A pseudo-R-squared corresponding to unity minus posterior variance over dependent variance. \cr
 \code{deviance.bma}: returns the \code{\link{deviance}} of a \code{bma} model as returned from \code{\link{bms}}. \cr
 \code{deviance.zlm}: returns the \code{\link{deviance}} of a \code{\link{zlm}} model.
}

\author{ Martin Feldkircher and Stefan Zeugner }
\seealso{ \code{\link{bms}} for creating \code{bma} objects and priors, \code{\link{zlm}} object.

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
 data(datafls)
  
 mm=bms(datafls[,1:10])
 deviance(mm)/nrow(datafls) # is equivalent to
 post.var(mm)
 
 post.pr2(mm) # is equivalent to
 1 - post.var(mm) / ( var(datafls[,1])*(1-1/nrow(datafls)) )
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
