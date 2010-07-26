\name{print.topmod}
\alias{print.topmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Printing topmod Objects }
\description{
  Print method for objects of class 'topmod', typically the best models stored in a 'bma' object
}
\usage{
 \method{print}{topmod}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class 'topmod' - see \code{\link{topmod}} }
  \item{\dots}{ additional arguments passed to \code{link{print}} }
}
\details{
  See \code{\link{pmp.bma}} for an explanation of likelihood vs. MCMC frequency concepts
}
\value{
  if \code{x} contains more than one model, then the function returns a 2-column matrix:
  \item{Row Names}{show the model binaries in hexcode (cf. \code{\link{topmodels.bma}})}
  \item{Column 'Marg.Log.Lik'}{shows the marginal log-likelihoods of the models in \code{x}}
  \item{Column 'MCMC Freq'}{shows the MCMC frequencies of the models in \code{x}}
    
  if \code{x} contains only one model, then more detailed information is shown for this model:
  \item{first line}{'Model Index' provides the model binary in hexcode, 'Marg.Log.Lik' its marginal log likelhood, 'Sampled Freq.' how often it was accepted (function \code{ncount()} in \code{\link{topmod}})}
  \item{Estimates}{first column: covariate indices included in the model, second column: posterior expected value of the coefficients, third column: their posterior standard deviations (excluded if no coefficients were stored in the topmod object - cf. argument \code{bbeta} in \code{\link{topmod}}) }
  \item{Included Covariates}{the model binary}
  \item{Additional Statistics}{any custom additional statistics saved with the model}
  
}
\author{ Martin Feldkircher and Stefan Zeugner }
\seealso{ \code{\link{topmod}} for creating topmod objects, \code{\link{bms}} for their typical use, \code{\link{pmp.bma}} for comparing posterior model probabilities 

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
# do some small-scale BMA for demonstration
data(datafls)
mm=bms(datafls[,1:10],nmodel=20)

#print info on the best 20 models
print(mm$topmod)
print(mm$topmod,digits=10)

#equivalent:
cbind(mm$topmod$lik(),mm$topmod$ncount())



#now print info only for the second-best model:
print(mm$topmod[2])

#compare 'Included Covariates' to:
topmodels.bma(mm[2])

#and to
as.vector(mm$topmod[2]$bool_binary())


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}