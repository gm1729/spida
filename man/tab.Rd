\name{tab}
\Rdversion{1.1}
\alias{tab}
\alias{tab.formula}
\alias{tab.data.frame}
\alias{tab.default}
\alias{Tab}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Table of frequencies or relative frequencies bordered with totals and including NAs
}
\description{
Generates a table of frequencies or relative frequencies or relative percentages
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{

tab( df, ~ a + b )

tab( df, ~ a + b , pct = 1)

tab( ~ a + b + c , df, pr = c(1,2))

tab( a, b, c )

tab( a, b, c, weights = w)

tab( w ~ a + b + c, df)

Tab( a, b, c, Weights = w, pct = 1, keep = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
as with \code{table}, one or more objects which can be interpreted as factors (including character strings),
or a list (or data frame) whose components can be so interpreted.
}
\item{data}{
  a data frame in which formula are interpreted
}
\item{fmla}{
   a formula whose right-hand side names the variables to be used for tabulation. The optional left-hand side specifies a variable to be used for weights.
}
\item{useNA}{
 whether to include NA levels. The default is "ifany". Can also be set to "no" or "always".
}
\item{pct}{
margins to be scaled to sum to 100
}
\item{pr}{
margins to be scaled to sum to 1
}
\item{total.margins}{
if FALSE, generate table without margins
}
\item{weights}{
instead of generating a frequency table, generate a table with the sum of the weights
}
\item{keep}{
names of margins to keep with 'Tab', default = "All". To drop all margins, use default = "".}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class 'table' of dimension equal to the number of variables, with optional margins showing totals. Elements of the matrix can be frequencies, relative frequencies, percentages or sums of weights. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Georges Monette
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
titanic <- as.data.frame(Titanic)
head(titanic)
tab(titanic, ~ Sex + Age, weights = Freq) 

tab(titanic, ~ Sex + Survived + Age, weights = Freq)

round(tab(titanic, ~ Sex + Survived + Age, 
    pct = c(1,3), weights = Freq),2)

round(Tab(titanic, ~ Sex + Survived + Age, 
    pct = c(1,3), weights = Freq),2)

round(Tab(titanic, ~ Sex + Survived + Age, 
    pct = c(1,3), weights = Freq, keep = ""),2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
