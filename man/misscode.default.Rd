\name{misscode.default}
\alias{misscode.default}
\alias{misscode.data.frame}
\alias{misscode.factor}
\alias{misscode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Turns NAs into a value below the range of non-missing data for plotting
%%  ~~function to do ... ~~
}
\description{
\code{misscode} turns NAs of numerical variables into a value slightly less that non-missing values. Missing values for factors are made into a non-missing level. When applied to a data frame, a new variable \code{.nmiss} is added to the data frame indicating the number of variables with missing data in each row of the data frame.
}
\usage{
misscode.default(x, ..., offset = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{offset}{
%%     ~~Describe \code{offset} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
function (x, ..., offset = 0.1) 
{
    rr <- range(x, na.rm = TRUE)
    vmiss <- min(x, na.rm = TRUE) - offset * diff(rr)
    nas <- is.na(x)
    x[nas] <- vmiss
    attr(x, "nas") <- nas
    x
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
