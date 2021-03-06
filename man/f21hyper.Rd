\name{f21hyper}
\alias{f21hyper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian Hypergeometric Function F(a,b,c,z)}
\description{
  Computes the value of a Gaussian hypergeometric function \eqn{ F(a,b,c,z) } for \eqn{-1 \leq z \leq 1} and \eqn{a,b,c \geq 0} 
}
\usage{
f21hyper(a, b, c, z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{ The parameter \code{a} of the Gaussian hypergeometric function, must be a positive scalar here }
  \item{b}{ The parameter \code{b} of the Gaussian hypergeometric function, must be a positive scalar here }
  \item{c}{ The parameter \code{c} of the Gaussian hypergeometric function, must be a positive scalar here  }
  \item{z}{ The parameter \code{z} of the Gaussian hypergeometric function, must be between -1 and 1 here }
}
\details{
  The function \code{f21hyper} complements the analysis of the 'hyper-g prior' introduced by Liang et al. (2008).\cr
  For parameter values, compare cf. \url{http://en.wikipedia.org/wiki/Hypergeometric_function#The_series_2F1}.
}
\value{
  The value of the Gaussian hypergeometric function \eqn{ F(a,b,c,z) }
}
\references{ Liang F., Paulo R., Molina G., Clyde M., Berger J.(2008): Mixtures of g-priors for Bayesian variable selection. J. Am. Statist. Assoc. 103, p. 410-423 

\url{http://en.wikipedia.org/wiki/Hypergeometric_function#The_series_2F1}
 }
\author{ Martin Feldkircher and Stefan Zeugner }
\note{ This function is a simple wrapper function of sped-up code that is intended for sporadic application by the user; it is neither efficient nor general; for a more general version cf. the package '\code{hypergeo}'}
\seealso{ package \code{hypergeo} for a more proficient implementation.

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
  
  f21hyper(30,1,20,.8) #returns about 165.8197
  
  f21hyper(30,10,20,0) #returns one
  
  f21hyper(10,15,20,-0.1) # returns about 0.4872972
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
