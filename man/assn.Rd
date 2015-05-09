\name{assn}
\Rdversion{1.1}
\alias{assn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assign in a pipeline
}
\description{
Same as \code{assign} with the first two arguments reversed and the default position referrring to the global environment so it can be used interactively in a pipeline
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
assn( value, name, pos = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{value}{
  to be assigned
}
  \item{x}{ name of the object to be assigned
}
  \item{pos}{ where the object is to be assigned. Default is \code{pos=0} which saves in the global environment
  }
%%     ~~Describe \code{B} here~~
  \item{...}{ other arguments are passed to \code{assign}
%%     ~~Describe \code{B} here~~
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
The formula is evaluated using \code{model.frame} and the result is used as the argument of \code{order} which, in turn is used to order the data frame. The ordered data frame is returned. 

For example, using the Prestige data set in package:car and the pipeline function in package:magrittr, we can sort the data frame and same the sorted version as "Prestige.ordered" with:

> require(car)
> require(magrittr)
> Prestige  %>% sortdf(~type+income)  %>% assn("Prestige.ordered")  %>%  dim 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
G. Monette
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

\dontrun{
require(magrittr)
require(car)
Prestige  %>% sortdf(~type+income)  %>% assn("Prestige.ordered")  %>%  dim 
}
}
