\name{osculant}
\alias{osculant.default}
\alias{osculant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Locus of osculation between two families of ellipses
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates points on the curve of osculation between the centers of two families of ellipses
}
\usage{
osculant(center1, shape1, center2, shape2, n = 100, range = c(0, 1), maxu = 100)

## Default S3 method
osculant(center1, shape1, center2, shape2, n = 100, range = c(0, 1), maxu = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{center1}{
%%     ~~Describe \code{center1} here~~
of first family of ellipses.
}
  \item{shape1}{
%%     ~~Describe \code{shape1} here~~
(i.e. 'variance matrix') of first family of ellipses.
}
  \item{center2}{
  of second family.
%%     ~~Describe \code{center2} here~~
}
  \item{shape2}{
%%     ~~Describe \code{shape2} here~~
of second family.
}
  \item{n}{
 n + 1 is the number of points to generate along the locus. \code{n = 1} generates the two centers, \code{n=0} generates to point on the first ellipse that lines on the locus of osculation provided the
 centre of the second family lines outside the first ellipse.
%%     ~~Describe \code{n} here~~
}
  \item{range}{
  of values of \code{u} to use to generate points. (See the algorithm in the code)
%%     ~~Describe \code{range} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The set of points on the locus is retuned as an \code{n+1} by \code{p} matrix where \code{p} is the
number of rows or columns of the shape matrices.
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
Georges Monette (georges@yorku.ca)
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{cell}}, \code{\link{dell}}, \code{\link{ellplus}},\code{\link{dellplus}},
}
\examples{
v1 <- 36*diag(2)
v2 <- .5 * (diag(2) - .4)
v2[2,2] <- 2
plot( 0:10,0:10,type = 'n')
lines( ell( c(2,2), v1))
lines( ell( c(4,4), v2), col = 'red')
osculant(  c(2,2), v1, c(4,4), v2, n = 3)
osculant(  c(2,2), v1, c(4,4), v2, n = 1)
lines( osculant( c(2,2), v1, c(4,4), v2), col = 'red')

lines( ell( c(8,8), v2), col = 'blue')
lines( osculant( c(2,2), v1, c(8,8), v2), col = 'blue')
points( osculant( c(2,2), v1, c(8,8), v2, n=1),pch = 16, col = 'blue')
points( osculant( c(2,2), v1, c(8,8), v2, n=0),pch = 16, col = 'blue')
points( osculant( c(8,8), v2, c(2,2), v1,  n=0),pch = 16, col = 'blue')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ellipse }
\keyword{ ellipse geometry }% __ONLY ONE__ keyword per line
