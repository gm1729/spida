\name{multresp}
\alias{multresp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create an indicator matrix from a multiple response variable
}
\description{
A multiple reponse variable is one whose values consist of a list of possible choices. Such a variable can be represented in a linear model through a matrix of indicator variables. 
}
\usage{
multresp(x, split = "[ \t,]+", sep = ".", root = deparse(substitute(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
a factor or character vector containing comma- or blank-separated lists of reponse choices for each subject. For example, each subject could have a list of symptoms: c('cough fever', 'cough bronchitis', 'headache wheeze', '', 'cough high_fever', NA). 
}
  \item{split}{
%%     ~~Describe \code{split} here~~
a regular expression matching the separator between fields. The default is \code{"[ \\t,]+"} which matches one more blank, tab or comma.  To use semicolons instead of commas as separators, use: \code{"[ \\t;]+"} .
}
  \item{sep}{
the separator used to generated the column names of the indicator matrix.  \code{sep} separates the root from the response values. Default: \code{"."}. For example, if \code{x} is the name of the variable, then the default column name for response value 'cough' is \code{x.cough}.
}
  \item{root}{
the root for column names. By default, the name of the variable.
the 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix of indicator variables. For, example, if 
\code{x <-c('cough fever', 'cough bronchitis', 'headache wheeze', '', 'cough high_fever', NA)}
\code{multresp(x)} produces:
\code{
     x.bronchitis x.cough x.fever x.headache x.high_fever x.wheeze
[1,]            0       1       1          0            0        0
[2,]            1       1       0          0            0        0
[3,]            0       0       0          1            0        1
[4,]            0       0       0          0            0        0
[5,]            0       1       0          0            1        0
[6,]            0       0       0          0            0        0

}
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
G. Monette (georges@yorku.ca)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x <-c('cough fever', 'cough bronchitis', 'headache wheeze', '', 'cough high_fever', NA)
multresp(x)

data <- data.frame(x = x, y = 1:6)
data <- with(data, cbind(data, multresp(x)))
data
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multiple response }
\keyword{ multiple response sets }% __ONLY ONE__ keyword per line
