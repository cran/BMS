\name{plotModelsize}
\alias{plotModelsize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Model Size Distribution}
\description{
  Plots posterior and prior model size distribution
}
\usage{
plotModelsize(bmao, exact = FALSE, ksubset = NULL, include.legend = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bmao}{ a 'bma' object (cf. \code{\link{bms}}) }
  \item{exact}{ if \code{TRUE}, then the posterior model distribution is based on the best models of \code{bmao} and their marginal likelihoods;\cr
                if \code{FALSE} (default) then the distribution is based on all encountered models and their MCMC frequencies (cf. 'Details' in \code{\link{coef.bma}}) }
  \item{ksubset}{ integer vector detailing for which model sizes the plot should be done }
  \item{include.legend}{ if TRUE, a small legend is included via the low-level command \code{\link{legend}} }
  \item{\dots}{ parameters passed on to \code{\link{matplot}} with sensible defaults  }
}
%\details{}
\value{
  As a default, \code{plotModelsize} plots the posterior model size distribution as a blue line, and the prior model distribution as a dashed red line.\cr
  In addition, it returns a list with the following elements:
  \item{mean}{The posterior expected value of model size}
  \item{var}{The variance of the posterior model size distribution}
  \item{dens}{A vector detailing the posterior model size distribution from model size \eqn{0} (the first element) to \eqn{K} (the last element)}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Martin Feldkircher and Stefan Zeugner }
\seealso{ See also \code{\link{bms}}, \code{\link{image.bma}}, \code{\link{density.bma}}, \code{\link{plotConv}}  

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
data(datafls)
mm=bms(datafls,burn=1500, iter=5000, nmodel=200,mprior="fixed",mprior.size=6)

#plot Nb.1 based on aggregate results
postdist= plotModelsize(mm)

#plot based only on 30 best models
plotModelsize(mm[1:30],exact=TRUE,include.legend=FALSE)

#plot based on all best models, but showing distribution only for model sizes 1 to 20
plotModelsize(mm,exact=TRUE,ksubset=1:20)

# create a plot similar to plot Nb. 1
plot(postdist$dens,type="l") 
lines(mm$mprior.info$mp.Kdist)


}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
