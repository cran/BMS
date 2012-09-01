\name{topmodels.bma}
\alias{topmodels.bma}
\title{Model Binaries and their Posterior model Probabilities}
%- Also NEED an '\alias' for EACH other topic documented here.
\description{
  Returns a matrix whose columns show which covariates were included in the best models in a 'bma' object. The last two columns detail posterior model probabilities.
}
\usage{
topmodels.bma(bmao)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bmao}{ an object of class 'bma' - see \code{\link{bms}} }
}
%\details{}
 
\value{
  Each column in the resulting matrix corresponds to one of the 'best' models in \code{bmao}: the first column for the best model, the second for the second-best model, etc.
  
  The model binaries have elements 1 if the regressor given by the row name was included in the respective models, and 0 otherwise.
  
  The second-last row shows the model's posterior model probability based on marginal likelihoods (i.e. its marginal likelihood over the sum of likelihoods of all best models)
  
  The last row shows the model's posterior model probability based on MCMC frequencies (i.e. how often the model was accepted vs sum of acceptance of all models)
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Martin Feldkircher and Stefan Zeugner }
\note{ Note that the column names are hexcode representations of the model binaries (e.g. "03" for \code{c(0,0,0,1,0,0)}) }
\seealso{ \code{\link{topmod}} for accessing model binaries, \code{\link{pmp.bma}} for posterior model probabilities 

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
data(datafls)
#sample with a limited data set for demonstration
mm=bms(datafls[,1:12],nmodel=20)

#show binaries for all
topmodels.bma(mm)

#show binaries for 2nd and 3rd best model, without the model probs
topmodels.bma(mm[2:3])[1:11,]

#access model binaries directly
mm$topmod$bool_binary()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

