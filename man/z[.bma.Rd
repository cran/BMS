\name{[.bma}
\alias{[.bma}
\alias{[.topmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract Parts of a topmod Object }
\description{
  returns a bma (or topmod) object whose 
}
\usage{
\method{[}{bma}(bmao, idx)

\method{[}{topmod}(tm, idx)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bmao}{ an object of class 'bma' }
  \item{idx}{ index of the best models to be retained }
  
  \item{tm}{ an object of class \code{\link{topmod}}' }
}
\details{
  Both funtions are equivalent: they serve to select a subset of the best models in a bma object.
}
\value{
  If \code{tm} is a topmod object (see \code{\link{topmod}}), then \code{tm[i]} returns a topmod object that only retains the models specified by \code{i}
  
  if \code{bmao} is an object of class bma (see \code{\link{bms}}), \code{bmao[i]} returns a bma object whose topmod object (i.e. \code{bmao$topmod} has beeen adjusted accordingly)
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
