\name{BMS-internal}
\title{Internal BMS objects}
\alias{.addup.enum}
\alias{.addup.enum.wotherstats}
\alias{.addup.mcmc}
\alias{.addup.mcmc.wotherstats}
\alias{.arglist}
\alias{.arglist.matchenvir}
\alias{.choose.gprior}
\alias{.choose.mprior}
\alias{.constr.intmat}
\alias{.construct.arglist}
\alias{.enum_fromindex}
\alias{.f21_4hyperg}
\alias{.f21simple}
\alias{.fls.samp}
\alias{.fls.samp.int}
\alias{.getpolycoefs}
\alias{.hexcode.binvec.convert}
\alias{.iterenum}
\alias{.iterenum.bone}
\alias{.iterenum.KgtN}
\alias{.lprob.constg.init}
\alias{.lprob.eblocal.init}
\alias{.lprob.hyperg.init}
\alias{.ols.terms2}
\alias{.post.beta.draws}
\alias{.post.calc}
\alias{.post.constant}
\alias{.post.estimates}
\alias{.post.topmod.bma}
\alias{.post.topmod.includes}
\alias{.rev.jump}
\alias{.rev.jump.int}
\alias{.starter}
\alias{.top10}
\alias{.topmod.as.bbetaT}
\description{Internal BMS objects.}
\details{These are not to be called by the user, but may come handy for programming tasks. 
There are some comments inside the functions to show how they work

The functions \code{.choose.mprior} and \code{.choose.gprior} may be interesting as they correspond to the elements \code{mprior.info} and \code{gprior.info} in a \code{link{bms}} result.

\code{.topmod.as.bbetaT} could be of use to a programmer interested in the \code{\link{topmod}} object.

\code{.ols.terms2} provides the basic OLS results on which this BMA routine is based.
The functions \code{.lprob.constg.init}, \code{.lprob.eblocal.init}, and \code{.lprob.hyperg.init} calculate g-prior-dependent posteriro statistics based on these OLS results.

}
\keyword{internal}