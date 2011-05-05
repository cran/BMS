\name{BMS-internal}
\title{Internal BMS objects}
\alias{.addup.enum}
\alias{.addup.enum.wotherstats}
\alias{.addup.mcmc}
\alias{.addup.mcmc.wotherstats}
\alias{.adjustdots}
\alias{.arglist}
\alias{.arglist.matchenvir}
\alias{.choose.gprior}
\alias{.choose.mprior}
\alias{.constr.intmat}
\alias{.construct.arglist}
\alias{.cor.topmod}
\alias{.enum_fromindex}
\alias{.f21_4hyperg}
\alias{.f21simple}
\alias{.fixedset.mprior}
\alias{.fixedset.sampler}
\alias{.fls.samp}
\alias{.getpolycoefs}
\alias{.hexcode.binvec.convert}
\alias{.iterenum}
\alias{.iterenum.bone}
\alias{.iterenum.KgtN}
\alias{.gprior.constg.init}
\alias{.gprior.eblocal.init}
\alias{.gprior.hyperg.init}
\alias{.mprior.customk.init}
\alias{.mprior.fixedt.init}
\alias{.mprior.pip.init}
\alias{.mprior.randomt.init}
\alias{.mprior.uniform.init}
\alias{.ols.terms2}
\alias{.post.beta.draws}
\alias{.post.calc}
\alias{.post.constant}
\alias{.post.estimates}
\alias{.post.topmod.bma}
\alias{.post.topmod.includes}
\alias{.quantile.density}
\alias{.rev.jump}
\alias{.starter}
\alias{.top10}
\alias{.topmod.as.bbetaT}
\description{Internal BMS objects.}
\details{These are not to be called by the user, but may come handy for programming tasks. 
The \code{BMS.RData} corresponding to the package (available at \url{http://bms.zeugner.eu} contains some comments inside the functions to show how they work

The functions \code{.mprior.fixedt.init} and \code{.gprior.constantg.init} are templates of funtions to create objects of class \code{\link{mprior-class}} and \code{\link{gprior-class}}, which may be used to customize \code{\link{bms}} and \code{\link{zlm}}.

\code{.topmod.as.bbetaT} could be of use to a programmer interested in the \code{\link{topmod}} object.

\code{.ols.terms2} provides the basic OLS results on which this BMA routine is based.

}
\keyword{internal}