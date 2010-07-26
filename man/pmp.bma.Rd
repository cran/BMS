\name{pmp.bma}
\alias{pmp.bma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Posterior Model Probabilities }
\description{
  Returns the posterior models probabilites based of the best models encountered by a 'bma' object
}
\usage{
pmp.bma(bmao)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bmao}{ A bma object (see argument \code{nmodel} in \code{\link{bms}}), alternatively an object of class 'topmod' }
}
\value{
 the result is a matrix, its column names describe the model binaries\cr
 There are two columns in the matrix:
  \item{PMP (Exact)}{posterior model probabilities based on the posterior likelihoods of the best models in \code{bmao}, relative to all 'best' models}
  \item{PMP (MCMC)}{posterior model probabilities of the best models in \code{bmao} based on their MCMC frequencies, relative to all models encountered by \code{bmao} - see 'Details' }
}
\details{
  A call to bms with an MCMC sampler (e.g. \code{bms(datafls,mcmc="bd",nmodel=100)} uses a Metropolis-Hastings algorithm to sample through the model space - and the frequency of how often models are drawn converges to the distribution of their posterior marginal likelihoods.
  While sampling, each 'bma' object stores the best models encountered by its sampling chain with their marginal likelihood and their MCMC frequencies.
  
  \code{pmp.bma} then allows for comparing the posterior model probabilities (PMPs) for the two different methods, similar to \code{\link{plotConv}}.
  It calculates the PMPs based on marginal likelihoods (first column) and the PMPs based on MCMC frequencies (second column) for the best x models stored in the bma object.
  
  The correlation of the two columns is an indicator of how well the MCMC sampler has converged to the actual PMP distribution - it is therefore also given in the output of \code{\link{summary.bma}}.
  
  The second column thus shows the PMPs of the best models relative to all models the call to \code{\link{bms}} has sampled through (therefore typically the second column adds up to less than one).
  However, the first column only relates to the likelihoods of the best models, therefore it adds up to 1. (the marginal likelihoods for the models not retained in the best models cannot be known).
  
  The second column is slightly different in case the \code{\link{bms}} argument \code{mcmc} was set to \code{mcmc="enumeration"}: In this case, the second column is also based on marginal likelihoods. The correlation between the two columns is therefore one.  
  
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Martin Feldkircher and Stefan Zeugner }
%\note{ ~~further notes~~ }

\seealso{ \code{\link{plotConv}} for plotting, \code{\link{bms}} for sampling bma objects 

Check \url{http://bms.zeugner.eu} for additional help.}

\examples{

## sample BMA for growth dataset, MCMC sampler
data(datafls)
mm=bms(datafls[,1:10],nmodel=20)

## mmodel likelihoods and MCMC frequencies of best 20 models
print(mm$topmod)

pmp.bma(mm)
#first column: posterior model prob based on model likelihoods,
#  relative to best models in 'mm'
#second column: posterior model prob based MCMC frequencies,
#  relative to all models encountered by 'mm'

#consequently, first column adds up to one
#second column shows how much of the sampled model space is
# contained in the best models
colSums(pmp.bma(mm))


#correlation betwwen the two shows how well the sampler converged
cor(pmp.bma(mm)[,1],pmp.bma(mm)[,2])

#is the same as given in summary.bma
summary(mm)["Corr PMP"]

#plot the two model probabilites
plotConv(mm)

#equivalent to the following chart
plot(pmp.bma(mm)[,2], type="s")
lines(pmp.bma(mm)[,1],col=2)


#moreover, note how the first column is constructed
liks=exp(mm$top$lik())
liks/sum(liks)
pmp.bma(mm)[,1] #these two are equivalent



#the example above does not converge well,
#too few iterations and best models
# this is already better, but also not good
mm=bms(datafls[,1:10],burn=2000,iter=5000,nmodel=200)


# in case the sampler has been 'enumeration' instead of MCMC,
# then both matrix columns are of course equivalent
mm=bms(datafls[,1:10],nmodel=512,mcmc="enumeration")
cor(pmp.bma(mm)[,1],pmp.bma(mm)[,2])
colSums(pmp.bma(mm))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
