\name{quantile.pred.density}
\alias{quantile.pred.density}
\alias{quantile.coef.density}
\alias{quantile.density}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract Quantiles from 'density' Objects }
\description{
  Quantiles for objects of class "density", "pred.density" or "coef.density"
}
\usage{
 \method{quantile}{pred.density}(x, probs = seq(.25,.75,.25), names = TRUE, ...)
 
 \method{quantile}{coef.density}(x, probs = seq(.25,.75,.25), names = TRUE, ...)
 
 \method{quantile}{density}(x, probs = seq(.25,.75,.25), names = TRUE, normalize = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a object of class \code{\link{pred.density}}, \code{coef.density}, \code{\link{density}}, or a list of densities. }
  \item{probs}{numeric vector of probabilities with values in [0,1] - elements very close to the boundaries return \code{Inf} or \code{-Inf}  }
  \item{names}{logical; if \code{TRUE}, the result has a \code{names} attribute, resp. a \code{rownames} and \code{colnames} attributes. Set to \code{FALSE} for speedup with many probs. }
  \item{normalize}{logical; if \code{TRUE} then the values in \code{x$y} are multiplied with a factor such that their integral is equal to one. }
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The methods \code{quantile.coef.density} and \code{quantile.pred.density} both apply \code{quantile.density} to densities nested with object of class \code{coef.density} or \code{pred.density}.\cr
  The function \code{quantile.density} applies generically to the built-in class \code{\link{density}} (as least for versions where there is no such method in the pre-configured packages).\cr
  Note that \code{quantile.density} relies on trapezoidal integration in order to compute the cumulative densities necessary for the calculation of quantiles. 
}

\value{
  If \code{x} is of class \code{density} (or a list with exactly one element), a vector with quantiles.\cr
  If \code{x} is a \code{\link{list}} of densities with more than one element (e.g. as resulting from \code{pred.density} or \code{coef.density}), then the output is a matrix of quantiles, with each matrix row corresponding to the respective density. 
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Stefan Zeugner }
\seealso{ 
\code{\link{quantile.default}} for a comparable function, \code{\link{pred.density}} and \code{\link{density.bma}} for the BMA-specific objects.

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
 data(datafls)
 mm = bms(datafls[1:70,], user.int=FALSE)
 
 #predict last two observations with preceding 70 obs:
 pmm = pred.density(mm, newdata=datafls[71:72,], plot=FALSE) 
 #'standard error' quantiles
 quantile(pmm, c(.05, .95))
 
 #Posterior density for Coefficient of "GDP60"
 cmm = density(mm, reg="GDP60", plot=FALSE) 
 quantile(cmm, probs=c(.05, .95))
 
 
 #application to generic density:
 dd1 = density(rnorm(1000))
 quantile(dd1)
 
 #application to list of densities:
 quantile.density( list(density(rnorm(1000)), density(rnorm(1000))) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
