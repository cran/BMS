\name{as.zlm}
\alias{as.zlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract a Model from a bma Object}
\description{Extracts a model out of a \code{bma} object's saved models and converts it to a \code{\link{zlm}} linear model}
\usage{
as.zlm(bmao, model = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bmao}{ A \code{bma} object, e.g. resulting from a call to \code{\link{bms}}}
  \item{model}{ The model index, in one of the following forms:\cr
An integer, denoting the rank of the model (1 for best, 2 for second-best, ...)\cr
A numeric or logical vector of length K describing which covariates are contained in the model\cr
A hexcode character describing which covariates are contained in the model
 }
}
\details{ A bma object stores several 'best' models it encounters (cf. argument \code{nmodel} in \code{\link{bms}}). \code{as.zlm} extracts a single model and converts it to an object of class \code{\link{zlm}}, which represents a linear model estimated under Zellner's g prior.}
\value{ a list of class \code{\link{zlm}}}
\author{
Stefan Zeugner
}
\seealso{
\code{\link{bms}} for creating \code{bma} objects, \code{\link{zlm}} for creating \code{zlm} objects, \code{\link{topmodels.bma}}  and \code{\link{pmp.bma}} for displaying the topmodels in a \code{bma} object

Check \url{http://bms.zeugner.eu} for additional help.
}
\examples{
data(datafls)

mm=bms(datafls[,1:6],mcmc="enumeration") # do a small BMA chain
topmodels.bma(mm)[,1:5] #display the best 5 models

m2a=as.zlm(mm,4) #extract the fourth best model
summary(m2a)

# extract the model only containing the 5th regressor
m2b=as.zlm(mm,c(0,0,0,0,1)) 

# extract the model only containing the 5th regressor in hexcode
print(bin2hex(c(0,0,0,0,1)))
m2c=as.zlm(mm,"01")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
