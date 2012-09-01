\name{pmpmodel}
\alias{pmpmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Posterior Model Probability for any Model }
\description{
  Returns the posterior model probability for any model based on bma results
}
\usage{
pmpmodel(bmao, model = numeric(0), exact = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bmao}{ A bma object as created by \code{\link{bms}}.}
  \item{model}{A model index - either variable names, or a logical with model binaries, or the model hexcode (cf. \code{\link{hex2bin}}, or a numeric with positions of the variables to be included.}
  \item{exact}{If \code{TRUE}, then the resulting PMP is based on analytical model likelihoods (works for any model). \cr
    If \code{FALSE}, the the resulting PMP is derived from MCMC frequencies (works only for the null and fullmodel, as well as for models contained in \code{bmao}'s topmod object.\cr
    If \code{bmao} is based on enumeration (cf. argument \code{mcmc} in \code{\link{bms}}, then \code{exact} does not matter. }
}
\value{
 A scalar with (an estimate of) the posterior model probability for \code{model}
}
\details{
  If the model as provided in \code{model} is the null or the full model, or is contained in \code{bmao}'s topmod object (cf. argument \code{nmodel} in \code{\link{bms}}), \cr
  then the result is the same as in \code{\link{pmp.bma}}.\cr
  If not and \code{exact=TRUE}, then \code{pmpmodel} estimates the model based on comparing its marginal likelihood (times model prior) to the likelihoods in the \code{topmod} object and multiplying by their sum of PMP according to MCMC frequencies,  
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Stefan Zeugner }
%\note{ ~~further notes~~ }

\seealso{ \code{\link{pmp.bma}} and \code{\link{topmodels.bma}} for similar functions 

Check \url{http://bms.zeugner.eu} for additional help.}

\examples{

## sample BMA for growth dataset, enumeration sampler
data(datafls)
mm=bms(datafls[,1:10],nmodel=5)

#show the best 5 models:
pmp.bma(mm)
#first column: posterior model prob based on model likelihoods,
#second column: posterior model prob based MCMC frequencies,

### Different ways to get the same result: #########

#PMP of 2nd-best model (hex-code representation)
pmpmodel(mm,"00c")

#PMP of 2nd-best model (binary representation)
incls=as.logical(beta.draws.bma(mm)[,2])
pmpmodel(mm,incls)

#PMP of 2nd-best model (via variable names)
#names of regressors in model "00c": 
names(datafls[,2:10])[incls]
pmpmodel(mm,c("SubSahara", "LatAmerica"))

#PMP of 2nd-best model (via positions)
pmpmodel(mm,c(6,7))

####PMP of another model #########
pmpmodel(mm,1:5)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
