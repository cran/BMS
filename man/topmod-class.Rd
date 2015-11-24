\name{topmod-class}
\docType{class}
\alias{topmod-class}

\title{Class "topmod"}
\description{ An updateable list keeping the best x models it encounters in any kind of model iteration }
\section{Objects from the Class}{
Objects can be created by calls to \code{\link{topmod}}, or indirectly by calls to \code{\link{bms}}.\cr

A 'topmod' object (as created by \code{topmod}) holds three basic vectors: \code{lik} (for the (log) likelihood of models or similar), \code{bool()} for a hexcode presentation of the model binaries (cf. \code{\link{bin2hex}})  and ncount() for the times the models have been drawn.\cr
All these vectors are sorted descendantly by \code{lik}, and are of the same length. The maximum length is limited by the argument \code{nbmodels}.
  
If \code{tmo} is a topmod object, then a call to \code{tmo$addmodel} (e.g. \code{tmo$addmodel(mylik=4,vec01=c(T,F,F,T))} updates the object \code{tmo} by a model represented by \code{vec01} (here the one including the first and fourth regressor) and the marginal (log) likelihood \code{lik} (here: 4).\cr
If this model is already part of \code{tmo}, then its respective \code{ncount} entry is incremented by one; else it is inserted into a position according to the ranking of \code{lik}.\cr
In addition, there is the possibility to save (the first moments of) coefficients of a model (\code{betas}) and their second moments (\code{betas2}), as well as an arbitrary vector of statistics per model (\code{fixed_vector}).\cr
	 
}
\section{Slots}{
  \describe{
    \item{\code{.S3Class}:}{Object of class \code{"list"}, elements are: 
     \describe{
    \item{\code{addmodel}:}{function that adjusts the list of models in the 'topmod' object (see Details). \code{mylik} is the basic selection criterion (usually log likelihood), \code{vec01} is the model binary (logical or numeric) indicating which regressors are included - cf. \code{\link{topmod}} }
    \item{\code{lik}:}{the function \code{lik()} returns a numeric vector of the best models (log) likelihoods, in decreasing order}
    \item{\code{bool}:}{the function \code{bool()} returns a character vector of hexmode expressions for the model binaries (cf. \code{\link{bin2hex}}), sorted by \code{lik()}  }
    \item{\code{ncount}:}{the function \code{ncount()} returns a numeric vector of MCMC frequencies for the best models (i.e. how often the respective model was introduced by \code{addmodel})}
    \item{\code{nbmodels}:}{the function \code{nbmodels()} returns the argument \code{nbmodel} to function \code{\link{topmod}} }
    \item{\code{nregs}:}{the function \code{nregs()} returns the argument \code{nmaxregressors} to \code{\link{bms}} }
    \item{\code{betas_raw}:}{ the function \code{betas_raw()}  returns a vector containing the coefficients in \code{betas} (see below) without the zero entries }
    \item{\code{betas2_raw}:}{ the function \code{betas2_raw()}  returns a vector containing the coefficient second moments in \code{betas2} (see below) without the zero entries }
    \item{\code{kvec_raw}:}{ the function \code{kvec_raw()} returns a vector with model sizes (integers) for the models denoted in \code{bool} }
    \item{\code{bool_binary}:}{the function \code{bool_binary()} returns a matrix whose columns present the models conforming to \code{lik()} in binary form }
    \item{\code{betas}:}{ the function \code{betas()} returns a matrix whose columns are the cofficents conforming to \code{bool_binary()} (Note that these include zero coefficents due to non-inclusion of covariates); Note: may be an empty matrix }
    \item{\code{betas2}:}{ the function \code{betas2()} returns a matrix similar to \code{betas()}, but with the coeffficents second moments (Note: can be empty)  }
    \item{\code{fixed_vector}:}{the function \code{fixed_vector()} returns a matrix whose columns bear the \code{fixed_vector} statistics conforming to \code{lik()} (see Details); Note: if \code{lengthfixedvec=0} in \link{topmod} this returns an empty matrix }
    }}
  }
}
\section{Methods}{
 \code{\link{print.topmod}}
}
\references{ \url{http://bms.zeugner.eu}  }
\author{Martin Feldkircher and Stefan Zeugner}
\seealso{
	\code{\link{topmod}} to create \code{topmod} objects and a more detailed description,\cr
    \code{\link{[.topmod}} for subselections, \code{\link{is.topmod}} to test for this class
}
\examples{
  tm= topmod(2,4,TRUE,0) #should keep a  maximum two models
  tm$addmodel(-2.3,c(1,1,1,1),1:4,5:8) #update with some model
  tm$addmodel(-2.2,c(0,1,1,1),1:3,5:7) #add another model
  tm$addmodel(-2.2,c(0,1,1,1),1:3,5:7) #add it again -> adjust ncount
  tm$addmodel(-2.5,c(1,0,0,1),1:2,5:6) #add another model
  
  #read out
  tm$lik()
  tm$ncount()
  tm$bool_binary()
  tm$betas()
}
\keyword{classes}
