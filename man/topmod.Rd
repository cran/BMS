% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_outer.R
\name{topmod}
\alias{topmod}
\alias{is.topmod}
\title{Topmodel Object}
\usage{
topmod(
  nbmodels,
  nmaxregressors = NA,
  bbeta = FALSE,
  lengthfixedvec = 0,
  liks = numeric(0),
  ncounts = numeric(0),
  modelbinaries = matrix(0, 0, 0),
  betas = matrix(0, 0, 0),
  betas2 = matrix(0, 0, 0),
  fixed_vector = matrix(0, 0, 0)
)
}
\arguments{
\item{nbmodels}{The maximum number of models to be retained by the topmod
object}

\item{nmaxregressors}{The maximum number of covariates the models in the
topmod object are allowed to have}

\item{bbeta}{if \code{bbeta=TRUE}, then first and second moments of model
coefficients are stored in addition to basic model statistics (Note: if
\code{bbeta<0} then only the first moments are saved)}

\item{lengthfixedvec}{The length of an optional fixed vector adhering to
each model (for instance R-squared, etc). If \code{lengthfixedvec=0} then no
additonal fixed vector will be stored.}

\item{liks}{optional vector of log-likelihoods to initialize topmod object
with (length must be \code{<=nbmodels}) - see example below}

\item{ncounts}{optional vector of MCMC frequencies to initialize topmod
object with (same length as \code{liks}) - see example below}

\item{modelbinaries}{optional matrix whose columns detail model binaries to
initialize topmod object with (same nb columns as \code{liks}, nb rows as
\code{nmaxregressors}) - see example below}

\item{betas}{optional matrix whose columns are coefficients to initialize
topmod object with (same dimensions as \code{modelbinaries}) - see example
below}

\item{betas2}{optional matrix whose columns are coefficients' second moments
to initialize topmod object with (same dimensions as \code{modelbinaries}) -
see example below}

\item{fixed_vector}{optional matrix whose columns are a fixed vector
initialize topmod object with (same \code{ncol} as \code{modelbinaries}) -
see example below}
}
\value{
a call to \code{topmod} returns a list of class "topmod" with the
following elements:
\item{addmodel(mylik,vec01,vbeta=numeric(0),vbeta2=numeric(0),fixedvec=numeric(0))}{function
that adjusts the list of models in the 'topmod' object (see Details).
\code{mylik} is the basic selection criterion (usually log likelihood),
\code{vec01} is the model binary (logical or numeric) indicating which
regressors are included.\cr \code{vbeta} is a vector of length equal to
\code{sum(vec01)}, contianing only the non-zero coefficients (only accounted
for if \code{bbeta!=FALSE}). \code{vbeta2} is a similar vector of second
moments etc. (only accounted for if \code{bbeta=TRUE}); \code{fixedvec} is
an arbitrary vector of length \code{lengthfixedvec} (see above)}
\item{lik()}{A numeric vector of the best models (log) likelihoods, in
decreasing order} \item{bool()}{A character vector of hexmode expressions
for the model binaries (cf. \code{\link{bin2hex}}), sorted by \code{lik()} }
\item{ncount()}{A numeric vector of MCMC frequencies for the best models
(i.e. how often the respective model was introduced by \code{addmodel})}
\item{nbmodels}{Returns the argument \code{nbmodel}} \item{nregs}{Returns
the argument \code{nmaxregressors}} \item{bool_binary()}{Returns a matrix
whose columns present the models conforming to \code{lik()} in binary form}
\item{betas()}{a matrix whose columns are the coefficients conforming to
\code{bool_binary()} (Note that these include zero coefficients due to
non-inclusion of covariates); Note: if \code{bbeta=FALSE} this returns an
empty matrix} \item{betas2()}{similar to \code{betas} , for the second
moments of coefficients Note: if \code{bbeta<=0}, this returns an empty
matrix} \item{fixed_vector()}{The columns of this matrix return the
\code{fixed_vector} statistics conforming to \code{lik()} (see Details);
Note: if \code{lengthfixedvec=0} this returns an empty matrix}
}
\description{
Create or use an updateable list keeping the best x models it encounters
(for advanced users)
}
\details{
A 'topmod' object (as created by \code{topmod}) holds three basic vectors:
\code{lik} (for the (log) likelihood of models or similar), \code{bool()}
for a hexcode presentation of the model binaries (cf. \code{\link{bin2hex}})
and ncount() for the times the models have been drawn.\cr All these vectors
are sorted descendantly by \code{lik}, and are of the same length. The
maximum length is limited by the argument \code{nbmodels}.

If \code{tmo} is a topmod object, then a call to \code{tmo$addmodel} (e.g.
\code{tmo$addmodel(mylik=4,vec01=c(T,F,F,T))} updates the object \code{tmo}
by a model represented by \code{vec01} (here the one including the first and
fourth regressor) and the marginal (log) likelihood \code{lik} (here: 4).

If this model is already part of \code{tmo}, then its respective
\code{ncount} entry is incremented by one; else it is inserted into a
position according to the ranking of \code{lik}.

In addition, there is the possibility to save (the first moments of)
coefficients of a model (\code{betas}) and their second moments
(\code{betas2}), as well as an arbitrary vector of statistics per model
(\code{fixed_vector}).

\code{is.topmod} returns \code{TRUE} if the argument is of class 'topmod'
}
\note{
\code{topmod} is rather intended as a building block for programming;
it has no direct application for a user of the BMS package.
}
\examples{

#standard use
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
  
  is.topmod(tm)
  
  #extract a topmod oobject only containing the second best model
  tm2=tm[2]
  
  
  
  #advanced: should return the same result as
  #initialize
  tm2= topmod(2,4,TRUE,0, liks = c(-2.2,-2.3), ncounts = c(2,1), 
          modelbinaries = cbind(c(0,1,1,1),c(1,1,1,1)), betas = cbind(0:3,1:4), 
          betas2 = cbind(c(0,5:7),5:8)) 

  #update 
  tm$addmodel(-2.5,c(1,0,0,1),1:2,5:6) #add another model
  
  #read out
  tm$lik()
  tm$ncount()
  tm$bool_binary()
  tm$betas()

}
\seealso{
the object resulting from \code{\link{bms}} includes an element of
class 'topmod'

Check \url{http://bms.zeugner.eu} for additional help.
}
\keyword{models}
