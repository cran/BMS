\name{plotComp}
\alias{plotComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compare Two or More bma Objects }
\description{
  Plots a comparison of posterior inclusion probabilites, coefficients or their standard deviation between various bma objects
}
\usage{
plotComp(..., varNr = NULL, comp = "PIP", exact = FALSE, 
   include.legend = TRUE, add.grid = TRUE, do.par = TRUE, cex.xaxis = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ one or more objects of class 'bma' to be compared. \code{plotComp} passes on any other parameters in \code{\dots} to \code{\link{matplot}}. }
  \item{varNr}{ optionally, covariate indices to be included in the plot, can be either integer vector or character vector - see examples }
  \item{comp}{ a character denoting what should be compared: \code{comp="PIP"} (default) for posterior inclusion probabilities, \code{comp="Post Mean"} for coefficients, \code{comp="Post SD"} for their standard deviations, \code{comp="Std Mean"} or standardized coefficients, or \code{comp="Std SD"} for standardized standard deviations }
  \item{exact}{ if \code{FALSE}, the statistics to be compared are based on aggregate bma statistics, if \code{TRUE}, they are based solely on the best models retained in the bma objects}
  \item{include.legend}{ whether to include a default legend in the plot (custom legends can be added with the command \code{\link{legend}}) }
  \item{add.grid}{ whether to add a \code{\link{grid}} to the plot }
  \item{do.par}{ whether to adjust \code{par("mar")} in order to fit in the tick labels on the x-axis }
  \item{cex.xaxis}{ font size scaling parameter for the x-axis - cf. argument \code{cex.axis} in \code{\link{par}} }
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
plotComp(mm1,mm2,varNr=1:4, col=c("black","red"))

#plot only coefficients for covariates 'GDP60 ' and 'LifeExp'
plotComp(mm1,mm2,varNr=c("GDP60", "LifeExp"),comp="Post Mean")



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
