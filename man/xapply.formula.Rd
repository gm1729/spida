\name{xapply.formula}
\alias{xapply.formula}
\alias{xapply}
\alias{xpandlist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply a function over a ragged array and return a data frame
}
\description{
Splits the data into subsets, computes summary statistics for each, and returns the result in a convenient form.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xapply(

xapply.formula(formula, data, FUN, ..., subset, na.action = na.omit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
 a two-sided formula: the lhs identifies the variable(s) that are the first argument of FUN (if there is more than one as in \code{cbind( x, y)}), then \code{FUN}
 is applied to each in turn. Note that this behaviour is different from that of \code{capply}, where \code{FUN} is applied to the resulting data frame. The rhs identifies
 the variables to be used to split the data into subsets.
}
  \item{data}{
a data frame.
}
  \item{FUN}{
a function.
}
  \item{\dots}{
additional arguments to the function
}
  \item{subset}{
a condition to subset the data frame.
}
  \item{na.action}{
determines whether to omit or include NAs in the grouping factors. Use \code{na.include} so NAs will form a distinct level.
%%     ~~Describe \code{na.action} here~~
}
}
\details{
\code{xapply} works like \code{\link{aggregate}} except that the result is returned as a vector in a data frame with the levels of the \code{by} variable
replicated according to the length of the result.

The intention in writing \code{xapply} was to facilitate the creation of 'prediction data frames' extending \code{expand.grid} so that the values of
a variable can easily depend on the value of a factor. For example, predicting weight from height it might be desired to have a range of heights for each
sex that is consistent with the conditional distribution.

Suppose a data frame \code{hw} contains variables Sex, height and weight. Instead of
\preformatted{
   > fit <- lm ( weight ~ height * Sex, hw)
   > pred <- expand.grid( Sex = levels(hw$Sex), height = quantile( hw$height, c(0,5,25,50,75,95,100)/100))
   > pred$weight <- predict( fit, pred)
   > xyplot( weight ~ height, pred, groups = Sex, type = 'b')
}
we can have:
\preformatted{
   > fit <- lm ( weight ~ height * Sex, hw)
   > pred <- expand.grid( Sex = levels(hw$Sex) )
   > pred <- merge( pred, xapply( height ~ Sex, hw, quantile, c(0,5,25,50,75,95,100)/100))
   > pred$weight <- predict( fit, pred)
   > xyplot( weight ~ height, pred, groups = Sex, type = 'b')
}

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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
