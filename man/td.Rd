\name{td}
\Rdversion{1.1}
\alias{td}
\alias{gd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set lattice parameters for multiple groups
}
\description{
Designed to easily sets lattice parameters for multiple groups. Setting parameters
before calling the lattice function allows parameters to be used consistently in the group key.
}
\usage{
td(new = FALSE, col = c("#0080ff", "#ff00ff", "darkgreen", "#ff0000", "orange", "#00ff00", "brown"), lty = 1:7, lwd = 1, pch = 1:7, cex = 0.8, font = 1, fill = "transparent", col.line = col, col.symbol = col, alpha = 1, alpha.line = alpha, alpha.symbol = alpha, len = 7, long = FALSE, record = FALSE, basecol = NULL, colsets = c("plot.symbol", "plot.line", "dot.symbol", "dot.line", "cloud.3d", "box.dot"), ...)

td(plot.symbol=list(col='red',pch=17))

td(plot.symbol=list(col='red',pch=17),plot.line=list(lwd=2))

gd(n=3)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{in 'gd' specifies the number of distinct colours to generate to distinguish groups. 'gd' uses 'latticeExtra' to set defaults for a ggplot2-like appearance. Default is n = 4
  }
  \item{new}{
If new = TRUE, open a new window, otherwise modify the existing active window, if there is one.
%%     ~~Describe \code{new} here~~
}
  \item{col}{}
  \item{lty}{}
  \item{lwd}{}
  \item{pch}{}
  \item{cex}{
%%     ~~Describe \code{col} here~~
vectors for graphical parameters for each level of the groups variable
}
  \item{font}{
%%     ~~Describe \code{font} here~~
}
  \item{fill}{
%%     ~~Describe \code{fill} here~~
}
  \item{col.line}{
%%     ~~Describe \code{col.line} here~~
line colour if 'groups' not given
}
  \item{col.symbol}{
%%     ~~Describe \code{col.symbol} here~~
symbol colour if 'groups' not given
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{alpha.line}{
%%     ~~Describe \code{alpha.line} here~~
}
  \item{alpha.symbol}{
graphical parameters for superpose.line and superpose.symbol
%%     ~~Describe \code{alpha.symbol} here~~
}
  \item{len}{
extend the length of parameters by recycling to lenght 'len'
%%     ~~Describe \code{len} here~~
}
  \item{long}{
  if TRUE generate a default combination of col, lty and pch with length 42.
%%     ~~Describe \code{long} here~~
}
  \item{record}{
  If TRUE, set history to 'recording'. Caution: this can use a lot of memory in .GlobalEnv.
  Consider adding selected graphs to memory with the 'Insert' key instead.
%%     ~~Describe \code{record} here~~
}
  \item{basecol}{
%%     ~~Describe \code{basecol} here~~
}
  \item{colsets}{
%%     ~~Describe \code{colsets} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
use to set any trellis parameter: e.g. \code{plot.symbol=list(cex=2,col='red')}. Particular useful for the cex, col, lty, lwd, alpha, pch parameters in plot.line and plot.symbol. 
}
}
\details{
'td' calls 'trellis.device' and sets graphical parameters for
'superpose.line' and 'superpose.symbol'. 'td' also initializes
a new trellis device with a white background if new = TRUE.
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

td( lty = 1:7)

td( plot.symbol = list(col = 'red', pch = 17))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
