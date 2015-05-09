\name{xqplot}
\Rdversion{1.1}
\alias{xqplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extended Quantile Plots
}
\description{
An easy way to see a dataset's variables at a glance. Shows uniform quantile plot for
numerical varibles and barcharts for factors. Quantile plots also show a horizontal line at the position of the mean and at mean plus or minus one standard deviation.
}
\usage{
xqplot(x, ptype = "quantile", labels = dimnames(x)[[2]], 
       ..., 
       mfrow = findmfrow(ncol(x)), 
       ask = prod(mfrow) < ncol(x) && dev.interac
       tive(), 
       mcex = 0.8, maxlab = 12, debug = F, 
       mar = F, maxvarnamelength = 20)
}
%- maybe alsx    : o 'usage' for other objects documented her c(5, 2, 3, 1), text.cex.factor = 1, left.labs =e.
\arguments{
  \item{x}{
  a data frame or list of variables to plot
%%     ~~Describe \code{x} here~~
}
  \item{ptype}{
     "quantile" (default) or "normal": kind of quantile to plot on x axis.
%%     ~~Describe \code{ptype} here~~
}
  \item{labels}{
  names for each plot
%%     ~~Describe \code{labels} here~~
}
  \item{\dots}{
  additional arguments passed to 'plot' command
%%     ~~Describe \code{\dots} here~~
}
  \item{mfrow}{
  number of rows and columns per page. If missing, an attempt
            is made to choose a reasonable number.
%%     ~~Describe \code{mfrow} here~~
}
  \item{ask}{
%%     ~~Describe \code{ask} here~~
}
  \item{mcex}{
  character expansion factor for marginal text
%%     ~~Describe \code{mcex} here~~
}
  \item{maxlab}{
  maximum number of categories to label in barcharts
%%     ~~Describe \code{maxlab} here~~
}
  \item{debug}{
  if TRUE, print additional information
%%     ~~Describe \code{debug} here~~
}
  \item{mar}{
  size of margins

%%     ~~Describe \code{mar} here~~
}
  \item{text.cex.factor}{
  character expansion factor for barchart labels
%%     ~~Describe \code{text.cex.factor} here~~
}
  \item{left.labs}{
  determines placement of barchart labels
%%     ~~Describe \code{left.labs} here~~
}
  \item{maxvarnamelength}{
  maximum length of variable name without splitting
            on two lines.
%%     ~~Describe \code{maxvarnamelength} here~~
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
Georges Monette <georges@yorku.ca>
}
\note{
Bugs:

     'mfrow' should take the total number of variables into account if they will
     fill more than one page so the last page is close to being full.

     The current version of the function could be made much simpler and
     more transparent. Some code is redundant.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(car)
xqplot(Prestige)
xqplot(Prestige,"n") # normal quantiles

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
