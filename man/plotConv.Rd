\name{plotConv}
\alias{plotConv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Convergence of BMA Sampler }
\description{
  Plots the posterior model probabilites based on 1) marginal likelihoods and 2) MCMC frequencies for the best models in a 'bma' object and details the sampler's convergence by their correlation
}
\usage{
plotConv(bmao, include.legend = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bmao}{ an object of class 'bma' - see \code{\link{bms}} }
  \item{include.legend}{ whether to include a legend in the plot }
  \item{\dots}{ other parameters for \code{\link{matplot}} }
}
\details{
    A call to bms with a MCMC sampler (e.g. \code{bms(datafls,mcmc="bd",nmodel=100)} uses a Metropolis-Hastings algorithm to sample through the model space: the frequency of how often models are drawn converges to the distribution of their posterior marginal likelihoods.\cr
    While sampling, each 'bma' object stores the best models encountered by its sampling chain with their marginal likelihood and their MCMC frequencies.\cr
  \code{plotConv} compares the MCMC frequencies to marginal likelihoods, and thus visualizes how well the sampler has converged.

}
\note{ \code{plotConv} is also used by \code{\link{plot.bma}} }
%\references{ ~put references to the literature/web site here ~ }
\author{Martin Feldkircher and Stefan Zeugner}
\seealso{ \code{\link{pmp.bma}} for posterior model probabilites based on the two concepts, \code{\link{bms}} for creating objects of class 'bma'

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{

data(datafls)
mm=bms(datafls[,1:12],user.int=FALSE)

plotConv(mm)

#is similar to
matplot(pmp.bma(mm),type="l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{aplot}

