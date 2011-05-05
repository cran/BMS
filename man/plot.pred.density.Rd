\name{plot.pred.density}
\alias{plot.pred.density}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Predictive Densities for bma Objects }
\description{
  Plots predictive densities for conditional forecasts of class 'pred.density' 
}
\usage{
## S3 method for class 'pred.density':
 \method{plot}{pred.density}(x, predict_index = NULL, addons = "eslz", realized.y = NULL, 
    addons.lwd = 1.5, ...) 
 }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{\link{pred.density}} }
  \item{predict_index}{An integer vector detailing which forecasted observations (corresponding to the argument \code{newdata} in \code{\link{pred.density}}) should be plotted.\cr
     Or the observation names to be plotted (as in \code{rownames(newdata)}). }
  \item{addons}{  character, defaulting to \code{"eslz"}. Specifies which additional information should be added to the plot via low-level commands (see 'Details' below).}
  \item{realized.y}{A vector with realized values of the dependent variables to be plotted in addition to the predictive density, must have its length conforming to \code{predict_index} (or \code{newdata})}
  \item{addons.lwd}{Line width to be used for the low-level plotting commands specified by \code{addons}. Cf. argument \code{lwd} in \code{\link{par}}.}
  \item{\dots}{arguments to be passed on to \code{\link{plot.density}}. }
}
\details{
The argument \code{addons} specifies what additional information should be added to the plot(s) via the low-level commands \code{\link{lines}} and \code{\link{legend}}:\cr
  "e" for the posterior expected value (EV) of the prediction,\cr
  "s" for 2 times its posterior standard deviation ('standard errors'),\cr
  "z" for a zero line, 
  "l" for including a \code{\link{legend}}\cr
  
  Any combination of these letters will give the desired result. Use \code{addons=""} for not using any of these.
  

}
%\references{ ~put references to the literature/web site here ~ }
\author{ Martin Feldkircher and Stefan Zeugner }
\seealso{ \code{\link{pred.density}} for constructing predictive densities, \code{\link{bms}} for creating \code{bma} objects, \code{\link{density.bma}} for plotting coefficient densities

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
 data(datafls)
 mm=bms(datafls,user.int=FALSE)
 
 #predictive density for two 'new' data points
 pd=pred.density(mm,newdata=datafls[1:2,]) 
 
 
 #plot the density for the second forecast observation
 plot(pd,2)  
 
 #plot the density with realized dep. variable, and no standard errors
 plot(pd, 1, realized.y=0,addons="lzeg") 
 
 
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{aplot}
