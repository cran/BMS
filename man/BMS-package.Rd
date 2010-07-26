\name{BMS-package}
\alias{BMS-package}
\alias{BMS}
\docType{package}
\title{
  Bayesian Model Sampling 0.2.4
}
\description{
This package enables Bayesian Model Averaging over the classical normal-conjugate model with many prior options and posterior statistics.
}
\details{
\tabular{ll}{
Package: \tab BMS\cr
Type: \tab Package\cr
Version: \tab 0.2.4\cr
Date: \tab 2010-07-26\cr
License: \tab Artistic 2.0\cr
}
The most important function is \code{\link{bms}} to perform bayesian model sampling.\cr
It basically offers to sample data according to different g-priors and model priors, and leaves the choice of different samplers (MCMC samplers, full or partial enumeration, and interaction samplers).\cr
The results provide analysis into models according to MCMC frequencies, and according to the posterior likelihood of the best \code{nmodel} models (cf. \code{\link{bms}}).

The functions \code{\link{coef.bma}} and \code{\link{summary.bma}} summarize the most important results.

The plotting functions \code{\link{plot.bma}}, \code{\link{image.bma}}, \code{\link{density.bma}}, and \code{\link{gdensity}} are the most important plotting functions (inter alia).

Moreover there are other functions for posterior results, such as \code{\link{beta.draws.bma}}, \code{\link{pmp.bma}}, and \code{\link{topmodels.bma}}, while \code{\link{c.bma}} helps to combine and parallelize sampling chains.

The function \code{\link{zlm}} estimates a Bayesian linear regression under Zellner's g prior, i.e. estimating a particular model without taking model uncertainty into account.

Finally, the small-scale functions \code{\link{f21hyper}}, \code{\link{hex2bin}} and \code{\link{fullmodel.ssq}} provide addidtional utilities.

Consider the functions \code{\link{topmod}} and as well as the internal functions \code{.choose.mprior} and \code{.choose.gprior} for more advanced programming tasks.
}
\author{
Martin Feldkircher and Stefan Zeugner

Maintainer: Stefan Zeugner <bms@zeugner.eu>
}
\references{
Feldkircher, M. and S. Zeugner (2009): Benchmark Priors Revisited: On Adaptive Shrinkage and the Supermodel Effect in Bayesian Model Averaging; IMF Working Paper 09-202
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
\seealso{
\url{http://bms.zeugner.eu}
}
\examples{
data(datafls)
mfls =bms(X.data=datafls,burn=1000,iter=9000,nmodel=100)
  info.bma(mfls)
  coef(mfls)
  coef(mfls,exact=TRUE,std.coefs=TRUE) 
  mfls[3]$topmod 
  image(mfls[1:20],FALSE) 
  plotModelsize(mfls,exact=TRUE) 
  density(mfls,"Spanish")
  
}
