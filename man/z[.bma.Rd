\name{[.bma}
\alias{[.bma}
\alias{[.topmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract Parts of a topmod Object }
\description{
  returns a bma (or topmod) object whose 
}
\usage{
\method{[}{bma}(x, i, ...)

\method{[}{topmod}(x, i, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class 'bma' or of class \code{\link{topmod}} }
  \item{i}{ index of the best models to be retained }
  \item{\ldots}{ arguments such as 'drop' will not have an effect }
}
\details{
  Both functions are equivalent: they serve to select a subset of the best models in a bma object.
}
\value{
  If \code{x} is a topmod object (see \code{\link{topmod}}), then \code{x[i]} returns a topmod object that only retains the models specified by \code{i}
  
  if \code{x} is an object of class bma (see \code{\link{bms}}), \code{x[i]} returns a bma object whose topmod object (i.e. \code{x$topmod} has been adjusted accordingly)
}
\author{ Martin Feldkircher and Stefan Zeugner }
\note{ Note that the method parameter \code{drop} is of no importance here. Moreover the corresponding assignment operator \code{\link{[<-}} is of no use.}
\seealso{ \code{\link{[}} for the default method, \code{\link{bms}} for crating bma objects 

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
#sample some BMA
 data(datafls)
 mm=bms(datafls[,1:11])
 
#plotConv for all models 
 plotConv(mm)

#plotConv for the best 20 models 
 plotConv(mm[1:20])

#model binary for the second best model 
topmodels.bma(mm[2])

#similar:
as.vector(mm$topmod[2]$bool_binary())

#exactly the same
as.vector(mm[2]$topmod$bool_binary())


#print stats for second best model
print(mm$topmod[2])




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
