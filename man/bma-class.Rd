\name{bma-class}
\docType{class}
\alias{bma-class}

\title{Class "bma"}
\description{ A list holding results from a BMA iteration chain }
\section{Objects from the Class}{
Objects can be created via calls to \code{\link{bms}}, but indirectly also via \code{\link{c.bma}}\cr
A \code{bma} object is a list whose elements hold information on input and output for a Bayesian Model Averaging iteration chain, such as from a call to \code{\link{bms}}:
}
\section{Slots}{
  \describe{
    \item{\code{info}:}{An object of class \code{"list"} holding aggregate statistics: \code{iter} is the number of iterations, \code{burn} the number of burn-ins. \cr
  The following have to be divided by \code{cumsumweights} to get posterior expected values: \code{inccount} are the posterior inclusion probabilities, \code{b1mo} and \code{b2mo} the first and second moment of coefficients, \code{add.otherstats} other statistics of interest (typically the moments of the shrinkage factor), \code{msize} is the post. expected model size, \code{k.vec} the posterior model size distribution, \code{pos.sign} the unconditional post. probability of positive coefficents, \code{corr.pmp} is the correlation between the best models' MCMC frequencies and their marg. likelihoods. \cr 
  \code{timed} is the time that was needed for MCMC sampling, \code{cons} is the posterior expected value of the constant. \code{K} and \code{N} are the maximum number of covariates and the sample size, respectively. }
    \item{\code{arguments}:}{An object of class \code{"list"} holding the evaluated function arguments provided to \code{\link{bms}}  }
    \item{\code{topmod}:}{An object of class \code{\linkS4class{topmod}} containing the best drawn models. see \code{\link{topmod}} for more details }
    \item{\code{start.pos}:}{the positions of the starting model. If bmao is a \code{bma} object this corresponds to covariates \code{bmao$reg.names[bmao$start.pos]}. If bmao is a chain that resulted from several starting models (cf. \code{\link{c.bma}}, then \code{start.pos} is a list detailing all of them.}
    \item{\code{gprior.info}:}{a list of class \code{\link{gprior-class}}, detailing information on the g-prior: \code{gtype} corresponds to argument \code{g} above, \code{is.constant} is FALSE if \code{gtype} is either "hyper" or "EBL", \code{return.g.stats} corresponds to argument \code{g.stats} above, \code{shrinkage.moments} contains the first and second moments of the shrinkage factor (only if \code{return.g.stats==TRUE}), \code{g} details the fixed g (if \code{is.constant==TRUE}), \code{hyper.parameter} corresponds to the hyper-g parameter \eqn{a} as in Liang et al. (2008).}
    \item{\code{mprior.info}:}{a list of class \code{\link{mprior-class}}, detailing information on the model prior: \code{origargs} lists the original arguments to \code{mprior} and \code{mprior.size} above; \code{mp.msize} denotes the prior mode size; \code{mp.Kdist} is a (K+1) vector with the prior model size distribution from 0 to K}
    \item{\code{X.data}:}{Object of class \code{"data.frame"} or class \code{"matrix"}, matrix: corresponds to argument \code{X.data} in \code{\link{bms}}, possibly cleaned for NAs }
    \item{\code{reg.names}:}{Vector of class \code{"character"}: the covariate names to be used for \code{X.data} in \code{\link{bms}} }
    \item{\code{bms.call}:}{Object of class \code{"call"}: the original call to the \code{\link{bms}} function }
  }
}
\section{Methods}{
 \code{\link{summary.bma}}, \code{\link{print.bma}}, \code{\link{coef.bma}}, \code{\link{density.bma}}, \code{\link{image.bma}}, \code{\link{plot.bma}}
}
\references{ \url{http://bms.zeugner.eu}  }
\author{Martin Feldkircher and Stefan Zeugner}
\seealso{
	\code{\link{bms}} for creating \code{bma} objects,\cr 
	or \code{\linkS4class{topmod}} for the topmod object
}
\examples{
 data(datafls)
 mm=bms(datafls)
 #show posterior model size
 print(mm$info$msize/mm$info$cumsumweights)
 #is the same number as in
 summary(mm)
 
}
\keyword{classes}
