\name{zlm-class}
\docType{class}
\alias{zlm-class}

\title{Class "zlm"}
\description{	A list holding output from the Bayesian Linar Model under Zellner's g prior akin to class 'lm' }
\section{Objects from the Class}{
Objects can be created via calls to \code{\link{zlm}}, but indirectly also via \code{\link{as.zlm}}.\cr
\code{\link{zlm}} estimates a Bayesian Linear Model under Zellner's g prior - its output is very similar to objects of class \code{\link{lm}} (cf. section 'Value')
}
\section{Slots}{
  \describe{
    \item{\code{coefficients}:}{A named vector of class \code{"numeric"}: holds the posterior expected values of 'regression' coefficients. The first element always holds the intercept}
    \item{\code{residuals}:}{Vector of class \code{"numeric"}: the residuals, that is the response variable minus fitted values }
    \item{\code{rank}:}{Scalar integer class \code{"numeric"}: the number of estimated parameters}
    \item{\code{fitted.values}:}{The (posterior expected values of the) fitted values}
    \item{\code{df.residual}:}{Scalar integer of class \code{"numeric"}: the residual degrees of freedom }
    \item{\code{call}:}{Object of class \code{"call"}:  the matched call to \code{\link{zlm}} that created the object}
    \item{\code{terms}:}{Object of class \code{"formula"}:  the \code{\link{terms}} object used}
    \item{\code{model}:}{Object of class \code{"data.frame"}: the model frame used }
    \item{\code{coef2moments}:}{Named vector of class \code{"numeric"}: coefficient posterior second moments }
    \item{\code{marg.lik}:}{Scalar integer of class \code{"numeric"}: the log marginal likelihood of the model }
    \item{\code{gprior.info}:}{An object of class \code{"list"} detailing information on the g-prior, cf. output value \code{gprior.info} in \code{\link{bms}} }
  }
}
\section{Methods}{
No methods defined with class "zlm" in the signature.
}
\references{ \url{http://bms.zeugner.eu}  }
\author{Martin Feldkircher and Stefan Zeugner}
\seealso{
	 \code{\link{zlm}} and \code{\link{as.zlm}} for creating \code{zlm} objects,\cr
     \code{\link{density.zlm}}, \code{\link{predict.zlm}} and \code{\link{summary.zlm}} for other posterior results 
}
%\examples{
%showClass("zlm")
%}
\keyword{classes}
