\name{tablemissing}
\alias{tablemissing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tabulate and visualize missing values
}
\description{
Tabulates missing values in variables and patterns of missing values among variables. It produces a missing value plot displaying different patterns of missing values among variables and missing value counts.
}
\usage{
tablemissing(x, sortby = "variable")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the dataset to be examined. It should be of class matrix or dataframe.
}
  \item{sortby}{
a logical value indicating whether the variables, patterns or both are sorted by the numbers of missing values. Possible values are "variable" (sort by variable name), "row" (sort by pattern counts), "column" (sort by variable counts) and "both" (sort by variable and pattern counts).
}
}
\details{
This function tabulates the missing values in variables and the patterns of missings. The input data should be of class matrix or dataframe. The function also displays a missing value plot. Each row represents a pattern of missings, where the variables with missing values are shown in yellow. The  barchart below shows the proportions of missing values by variable. The barchart on the right displays the pattern frequencies. This barchart uses the concept of censored zooming. The topmost bar, representing no missing values in any variable, is likely to be much larger than the others and will then not be drawn to scale, but bordered in red to show this. Ordering of variables and patterns of missings based on counts can be very useful for exploring the structure of missings in a dataset.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{finaltable }{A table of class dataframe containing number of missing values in variable and missing pattern}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Waqas Ahmed Malik
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(adult, package="vmv")
tablemissing(adult, sortby="variable")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ missing values }
\keyword{ plot}% __ONLY ONE__ keyword per line
