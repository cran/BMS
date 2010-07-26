\name{plotComp}
\alias{plotComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compare Two or More bma Objects }
\description{
  Plots a comparison of posterior inclusion probabilites, coefficients or their standard deviation between various bma objects
}
\usage{
plotComp(..., varNr = NULL, comp = "PIP", exact = FALSE, 
   cex.axis = 0.45, main = NULL, type = "p", lty = 1:5, lwd = 1.5, pch = NULL, 
   col = NULL, cex = NULL, bg = NA, xlab = "", ylab = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ two or more objects of class 'bma' }
  \item{varNr}{ optionally, covariate indices to be included in the plot, can be either integer vector or character vector - see examples }
  \item{comp}{ a character denoting what should be compared: \code{comp="PIP"} (default) for posterior inclusion probabilities, \code{comp="Post Mean"} for coefficients, \code{comp="Post SD"} for their standard deviations, \code{comp="Std Mean"} or standardized coefficients, or \code{comp="Std SD"} for standardized standard deviations }
  \item{exact}{ if \code{FALSE}, the statistics to be compared are based on aggregate bma statistics, if \code{TRUE}, they are based solely on the best models retained in the bma objects}
  \item{cex.axis}{ font size parameter for the axes - see \code{\link{par}} }
  \item{main}{ parameter for \code{\link{matplot}} - cf. \code{link{plot.default}} }
  \item{type}{ parameter for \code{\link{matplot}}, default \code{"p"} }
  \item{lty}{ parameter for \code{\link{matplot}}, default \code{1:5} }
  \item{lwd}{ parameter for \code{\link{matplot}}, default \code{1} }
  \item{pch}{ parameter for \code{\link{matplot}}, default \code{NULL} }
  \item{col}{ parameter for \code{\link{matplot}}, default special colors }
  \item{cex}{ parameter for \code{\link{matplot}}, default \code{NULL} }
  \item{bg}{ parameter for \code{\link{matplot}}, default \code{NA} }
  \item{xlab}{ parameter for \code{\link{matplot}}, default \code{""} }
  \item{ylab}{ parameter for \code{\link{matplot}}, default argument \code{comp} }

}

\author{ Martin Feldkircher and Stefan Zeugner }
\seealso{ \code{\link{coef.bma}} for the underlying function 

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
## sample two simple bma objects
data(datafls)
mm1=bms(datafls[,1:15])
mm2=bms(datafls[,1:15])

#compare PIPs
plotComp(mm1,mm2)

#compare standardized coefficeitns
plotComp(mm1,mm2,comp="Std Mean")

#...based on the lieklihoods of best models 
plotComp(mm1,mm2,comp="Std Mean",exact=TRUE)

#plot only PIPs for first four covariates
plotComp(mm1,mm2,varNr=1:4)

#plot only coefficients for covariates 'GDP60 ' and 'LifeExp'
plotComp(mm1,mm2,varNr=c("GDP60", "LifeExp"),comp="Post Mean")



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
