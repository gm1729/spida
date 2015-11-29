\name{dropone.lme}
\alias{dropone.lme}
\alias{dropone}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Refit a model dropping one observation or one cluster at a time
%%  ~~function to do ... ~~
}
\description{
Refits a model droping one observation or cluster at a time and returns a 
data frame with the sesulting fitted parameters.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dropone(fit, by, verbose = FALSE, ...)


dropone.lme(fit, by, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
The fitted object from which observations are dropped.
%%     ~~Describe \code{fit} here~~
}
  \item{by}{
an optional formula giving the variable to be used to drop observations. If missing observations are dropped one at at time.
%%     ~~Describe \code{by} here~~
}
\item{verbose}{
produce verbose output (default FALSE).
%%     ~~Describe \code{verbose} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame with cluster-level variables and the values of 'drop one' estimated parameters.
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
function (fit, by, verbose = FALSE, ...) 
{
    getSummary <- function(x) {
        re <- getVarCov(x)
        rn <- expand.grid(dimnames(re))
        rn <- array(apply(rn, 1, paste, collapse = "+"), dim = dim(re))
        re <- re[col(re) <= row(re)]
        names(re) <- rn[col(rn) <= row(rn)]
        ff <- fixef(x)
        names(ff) <- paste("b", names(ff), sep = ".")
        c(ff, dvcov = det(vcov(x)), re, sigma2 = fit$sigma^2)
    }
    data <- getData(fit)
    if (missing(by)) 
        by <- 1:nrow(data)
    if (inherits(by, "formula")) {
        by <- model.frame(by, data, na.action = na.include)
        by <- apply(by, 1, paste, collapse = "/")
    }
    by <- as.character(by)
    levs <- unique(by)
    data$.drop <- by
    data.ret <- up(data, ~.drop)
    names(levs) <- levs
    ret <- lapply(levs, function(x) {
        try(getSummary(update(fit, data = data[by != x, ])))
    })
    ret <- lapply(ret, function(x) if (inherits(x, "try-error")) 
        NA
    else x)
    ret <- do.call(rbind, ret)
    disp(dim(ret))
    ret <- as.data.frame(ret)
    ret$.drop <- levs
    disp(dim(ret))
    disp(dim(data.ret))
    ret <- merge(data.ret, ret, by = ".drop")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
