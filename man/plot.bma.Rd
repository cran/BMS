\name{plot.bma}
\alias{plot.bma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Posterior Model Size and Model Probabilities }
\description{
  Produces a combined plot: upper row shows prior and posterior model size distribution, lower row shows posterior model probabilities for the best models 
}
\usage{
 \method{plot}{bma}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class 'bma' }
  \item{\dots}{ additional arguments for \code{\link{matplot}} }
}
\value{
  combines the plotting functions \code{\link{plotModelsize}} and \code{\link{plotConv}}
}
\author{ Martin Feldkircher and Stefan Zeugner }
\note{ The upper plot shows the prior and posterior distribution of model sizes (\code{\link{plotModelsize}}).\cr
The lower plot is an indicator of how well the bma object has converged (\code{\link{plotConv}}). }
\seealso{ \code{\link{plotModelsize}} and \code{\link{plotConv}} 

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
data(datafls)
mm=bms(datafls,user.int=FALSE)

plot(mm)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
