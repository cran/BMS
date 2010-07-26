\name{bin2hex}
\alias{bin2hex}
\alias{hex2bin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Converting Binary Code to and from Hexadecimal Code }
\description{
  A simple-to-use function for converting a logical ('binary') vector into
  hex code and reverse.
}
\usage{
bin2hex(binvec)

hex2bin(hexcode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{binvec}{a logical vector (alternatively a vector coercible into logical)}
  \item{hexcode}{a single-element character denoting an integer in hexcode (admissible character: 0 to 9, ato f)}
}
\details{
  The argument is an integer in binary form (such as "101"), provided as a logical (\code{c(T,F,T)}) or numeric vector (\code{c(1,0,1)}).\cr 
  \code{bin2hex} then returns a character denoting this number in hexcode (in this case "5").
  
  The function \code{hex2bin} does the reverse operation, e.g. \code{hex2bin("5")} gives (\code{c(1,0,1)}).
  }
\value{
  \code{bin2hex} returns a single element character; \code{hex2bin} returns a numeric vector equivalent to a logical vector
}

\author{ Martin Feldkircher and Stefan Zeugner }
%\note{ Note that both functions preserve leading zeroes: for instance  \code{bin2hex(c(0,1))} returns a result different from \code{bin2hex(c(rep(0,8),1))}, and \code{hex2bin("0b")} differs from \code{hex2bin("b")} }
\seealso{ \code{\link{hex2bin}} for converting hexcode into binary vectors, \code{\link{format.hexmode}} for a related R function.

Check \url{http://bms.zeugner.eu} for additional help.}
\examples{
  bin2hex(c(TRUE,FALSE,TRUE,FALSE,TRUE,TRUE))
  bin2hex(c(1,0,1,0,1,1))
  hex2bin("b8a")
  bin2hex(hex2bin("b8a"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{arith}
