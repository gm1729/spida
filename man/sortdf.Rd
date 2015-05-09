\name{sortdf}
\Rdversion{1.1}
\alias{sortdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Order the rows of a data frame
}
\description{
Order the rows of a data frame using variables identified in a formula evaluated in the data frame. Convenient in a pipeline.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sordf( data, form )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  a data frame
%%     ~~Describe \code{A} here~~
}
  \item{form}{
  a formula identifying variables to be used in sorting
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

For example, using the Prestige data set in package:car and the pipeline function in package:magrittr, we can sort the data frame with:

> require(car)
> require(magrittr)
> Prestige  %>% sortdf(~type+income) -> Prestige. 
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
\dontrun{
require(car)
require(magrittr)
Prestige  %>% sortdf(~type+income)
}
}
