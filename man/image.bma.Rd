\name{image.bma}
\alias{image.bma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Signs of Best Models }
\description{
  Plots a grid with signs and inclusion of coefficients vs. posterior model probabilities for the best models in a 'bma' object: 
}
\usage{
 \method{image}{bma}(x, yprop2pip = FALSE, order.by.pip = TRUE, 
  do.par = TRUE, do.grid = TRUE, do.axis = TRUE, cex.axis = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a list of class bma (cf. \code{\link{bms}} for further details) }
  \item{yprop2pip}{ if \code{yprop2pip=TRUE} then the grid lines on the vertical axis are scaled according to the coefficients' inclusion probabilites.\cr
      If \code{yprop2pip=FALSE} (default) then the grid lines on the vertical axis are equidistant. 
      }
  \item{order.by.pip}{ with \code{order.by.pip=TRUE} (default), coefficients are sorted according to their posterior inclusion probabilites along the vertical axis. If \code{order.by.pip=FALSE} they are ordered as they were provided to \code{\link{bms}}.}
  \item{do.par}{ Defaults to \code{do.par=TRUE}, which adjusts \code{\link{par}()$mar} for optimal positioning. Set \code{do.par=FALSE} for customizing \code{par} yourself.}
  \item{do.grid}{ \code{do.grid=TRUE} (default) plots grid lines among the chart's boxes, akin to the low level command \code{\link{grid}}. \code{do.grid=FALSE} omits the grid lines.}
  \item{do.axis}{ \code{do.axis=TRUE} (default) plots axis tick marks and labels (cf. \code{\link{axis}}). \code{do.axis=FALSE} omits them. }
  \item{cex.axis}{ font size for the axes (cf. \code{\link{axis}}), defaults to 1  }
  \item{\dots}{ Parameters to be passed on to \code{\link{image.default}}.}
}
\details{
  Under default settings, blue corresponds to positive sign, red to a negative sign, white to non-inclusion.
}
%\value{
%  ~Describe the value returned
%}
%\references{ ~put references to the literature/web site here ~ }
\author{ Martin Feldkircher and Stefan Zeugner}
%\note{ ~~further notes~~ 
 %~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \link{coef.bma} for the coefficients in matrix form, \link{bms} for creating 'bma' objects. 

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
 data(datafls)
 
 model=bms(datafls,nmodel=200)
 
 #plot all models
 image(model,order.by.pip=FALSE)
 image(model,order.by.pip=TRUE,cex.axis=.8)
 
 #plot best 7 models, with other colors
 image(model[1:7],yprop2pip=TRUE,col=c("black","lightgrey"))
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
