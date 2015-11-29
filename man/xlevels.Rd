\name{xlevels}
\alias{xlevels}
\alias{dlevels}
\alias{Pmat}
\alias{Pmat_diffmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create 'derivatives' and 'means' of factors to generate, for example, pairwise differences or centres of existing factor for prediction
}
\description{
The functions \code{xlevels} and \code{dlevels} are primarily intended to create arguments for \code{expand.grid}
to create 'prediction' or 'effect data frames',
to generate wald tests and estimates of specific effects in models with interactions,
}
\usage{

xlevels(f, type = c("raw", "(mean)"), all = FALSE, sep = "")

dlevels(f, type = "pairwise", all = FALSE, sep = "")

Pmat(f, type = c("factor", "raw", "mean", "(mean)", "<mean>", "II", "cen", "cent", "center", "centre", "<centre>", "<center>", "(center)", "(centre)", "III", "pairwise", "<others.m>", "<others.c>", "diff", "diffmean", "diffc", "diffcen", "diffcentre", "diffcentre"), all = FALSE, sep = "")

Pmat_diffmat(w)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
a factor or otherwise a vector interpreted as the levels of a factor
%%     ~~Describe \code{f} here~~
}
  \item{type}{
one or more types of contrasts ('derivative') or means (convex combinations) of factor levels
}
  \item{all}{
for some values of \code{type}, indicates whether to use all contrasts. e.g. \code{type = "pairwise"} will produce pairwise comparisons
in both directions if \code{all == TRUE}
}
  \item{sep}{
can add spaces in constructed factor levels
}
  \item{w}{
weights used for each factor level in creating contrasts differentiating a factor level from others. Weights corresponding to
frequencies in the data frame result in effects corresponding to type II effects while equal weights correspond to type III effects for interacting
specific effects
%%     ~~Describe \code{w} here~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
