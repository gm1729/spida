\name{trim}
\Rdversion{1.1}
\alias{trim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trim trailing blanks
}
\description{
Generic function to trim leading and trailing blanks from character vectors and factors.
}
\usage{
trim(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a data frame, factor, character or numeric vector
%%     ~~Describe \code{x} here~~
}
}
\details{
The main application is in reading SPSS files that often have leading or trailing blanks in
character and factor values. These blanks are often inconsistent so that values will
appear to differ even though they are equal.  The trim function is called in \code{Read.spss} to
remove leading and trailing blanks from all factors.
}
\value{
A character vector or factor with leading and trailing blanks removed.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
trim in fun.R
    UseMethod("trim")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
