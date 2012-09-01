\name{variable.names.bma}
\alias{variable.names.bma}
\alias{model.frame.bma}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Variable names and design matrix }
\description{
  Simple utilities retrieving variable names and design matrix from a bma object
}
\usage{

 \method{variable.names}{bma}(object,...)
 \method{model.frame}{bma}(formula,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{bma} object (as produced by \code{\link{bms}}) }
  \item{formula}{ A \code{bma} object (as produced by \code{\link{bms}}) }
  \item{\ldots}{further arguments passed to or from other methods }
}
\details{
  All functions are \code{bma}-functions for the generic methods \code{\link{variable.names}}, \code{\link{deviance}}, and \code{\link{model.frame}}.
}

\author{ Martin Feldkircher and Stefan Zeugner }
\seealso{ \code{\link{bms}} for creating bma objects 

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
 data(datafls)
 bma_enum=bms(datafls[1:20,1:10])
 
 model.frame(bma_enum) # similar to 
 bma_enum$X.data
 
 variable.names(bma_enum)[-1] # is equivalent to
 bma_enum$reg.names
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
