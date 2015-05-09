\name{wald}
\Rdversion{1.1}
\alias{wald}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
General Linear Hypothesis for the 'fixed' portion of a model with Wald test
}
\description{
Tests a general linear hypothesis for the linear fixed portion of a model. The hypothesis
can be specified in a variety of ways such as a hypothesis matrix or a pattern that is used as
a regular expression to be matched with the names of coefficients of the model.
A number of tools are available to facilitate the generation of hypothesis matrices.
}
\usage{
wald(fit, Llist = "", clevel = 0.95, data = NULL, debug = FALSE, maxrows = 25,
full = FALSE, fixed = FALSE, invert = FALSE, method = "svd", help = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{ a model for which a \code{getFix} method exists.
%%     ~~Describe \code{fit} here~~
}
  \item{Llist}{ a hypothesis matrix or a pattern to be matched or a list of these
%%     ~~Describe \code{Llist} here~~
}
  \item{clevel}{ level for confidence intervals
%%     ~~Describe \code{clevel} here~~
}
  \item{data}{ used for 'data' attribute of value returned
%%     ~~Describe \code{data} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
}
  \item{maxrows}{ maximum number of rows of hypothesis matrix for which a full variance-covariance matrix is returned
%%     ~~Describe \code{maxrows} here~~
}
  \item{full}{ if TRUE, the hypothesis matrix is the model matrix for \code{fit} such that
  the estimated coefficients are the predicted values for the fixed portion of the model. This is designed
  to allow the calculation of standard errors for models for which the \code{predict} method does not provide them.
%%     ~~Describe \code{full} here~~
}
  \item{fixed}{ if \code{Llist} is a character to be used a regular expression, if \code{fixed} is TRUE
  \code{Llist} is interpreted literally, i.e. characters that have a special meaning in regular expressions
  are interpreted literally.
%%     ~~Describe \code{fixed} here~~
}
  \item{invert}{ if \code{Llist} is a character to be used a regular expression, \code{invert == TRUE}
  causes the matches to be inverted so that coefficients that do not match will be selected.
%%     ~~Describe \code{invert} here~~
}
  \item{method}{
     'svd' (current default) or 'qr' is the method used to find the full rank version of the hypothesis matrix.
     'svd' has correctly identified the rank of a large hypothesis matrix where 'qr' has failed.
}
  \item{help}{  obsolete
%%     ~~Describe \code{help} here~~
}
}
\details{
\verb{
General Linear Hypothesis with Wald test
      for lm, glm, lme, nlme and lmer objects.
      Can be extended to other objects (e.g.) 'glm' by writing 'getFix.glm'

      Usage:
      
         wald ( fit, L ) where L is a hypothesis matrix

         wald ( fit, 'pat' ) where 'pat' a regular expression (see ?regex)
           used to match names of coefficients of fixed effects.
           e.g. wald( fit, ':.*:') tests all 2nd and higher order interactions.

         wald ( fit, c(2,5,6)) to test 2nd, 5th and 6th coefficients.

         wald ( fit, list( hyp1= c(2,5,6), H2 = 'pat')) for more than one hypothesis
           matrix

      There are number of functions to help construct hypothesis matrices:

         Lform ( fit, list(expr1, expr2, ..., exprk), data = dframe)
           creates an L matrix by evaluating expr1, expr2, ..., exprk in the dframe
           environment to generate the k columns of the L matrix. 'dframe' is
           model.frame(fit) by default. See the example below to estimate a derivative.

      Creates an L matrix using formulas evaluated in 'data' for
      each column of the L matrix
      Example:

      library(car)
      fit <- lm( income ~ (education + I(education^2) )* type, Prestige)
      summary(fit)

      . . .
      Coefficients:
                               Estimate Std. Error t value Pr(>|t|)
       (Intercept)                891.3    23889.1   0.037  0.97032
       education                  210.0     5638.8   0.037  0.97037
       I(education^2)              38.3      328.3   0.117  0.90740
       typeprof                191523.2    63022.0   3.039  0.00312 **
       typewc                   25692.3    73888.0   0.348  0.72887
       education:typeprof      -28133.0    10236.0  -2.748  0.00725 **
       education:typewc         -4485.4    14007.9  -0.320  0.74956
       I(education^2):typeprof   1017.5      451.8   2.252  0.02679 *
       I(education^2):typewc      170.9      671.6   0.255  0.79967
      . . .

      # estimate the marginal value of education for each 
      # occupation in the data set

      L <- list( 'marginal value of education' =Lform( fit,
             form = list(0, 1, 2 * education, 0, 0, 
                          type == 'prof', 
                          type == 'wc',
                          2 * education * (type == 'prof'), 
                          2 * education * (type == 'wc')),
             data = Prestige))
      wald( fit, L )
      chat <- coef( wald( fit, L ), se = 2)
      xyplot( coef +coefp+coefm ~ education | type, cbind(Prestige,chat)[order(Prestige$education),],
              type = 'l')
      xyplot( chat~ education | type, Prestige)
      
      
      # This approach can be used to predict responses with a fitting method that has a 
      # 'model.matrix' method but does not have a 'predict' method or does not return 
      # estimated standard errors with its 'predict' method.
      
      datafit <- model.frame(fit)   # or just the data frame used for the fit
      ww <- wald(fit, model.matrix(fit))
      datafit <- cbind(datafit, coef(ww, se =2))
      # ...etc as above

      # To extend the 'wald' function to a new class of objects, one needs to write a 'getFix'
      # method to extract estimated coefficients, their estimated covariance matrix, and the
      # denominator degrees of freedom for each estimated coefficient. The getFix method for 
      # glm objects is:

      getFix.glm <- function(fit,...) {
             ss <- summary(fit)
             ret <- list()
             ret$fixed <- coef(fit)
             ret$vcov <- vcov(fit)
             ret$df <- rep(ss$df.residual, length(ret$fixed))
             ret
      }

      # and for 'mipo' objects in the packages 'mice':

      getFix.mipo <- function( fit, ...){ 
        #
        # pooled multiple imputation object in mice
        # 'wald' will use the minimal df for components with non-zero weights
        # -- this is probably too conservative and should be improved
        #
        ret <- list()
        ret$fixed <- fit$qbar
        ret$vcov <- fit$t
        ret$df <- fit$df
        ret
      }
}
}
\value{
An object of class \code{wald}, with the following components
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
\code{\link{Lform}},
\code{\link{xlevels}},
\code{\link{dlevels}},
\code{\link{Lall}},\code{\link{Lc}},\code{\link{Lequal}},
\code{\link{Ldiff}},\code{\link{Lmu}},\code{\link{Lmat}},\code{\link{Lrm}},
\code{\link{Leff}},
\code{\link{as.data.frame.wald}}. To extend to
new models see \code{\link{getFix}}. To generate
hypothesis matrices for general splines see \code{\link{gsp}} and \code{\link{sc}}.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
# The derivative with respect
# to ses could be evaluated at each point in the following:

data(hs)
require( nlme )
fit <- lme( mathach ~ (ses + I(ses^2)) * Sex, hs, random = ~ 1 + ses | school)

wald( fit, 'Sex')  # sig. overall effect of Sex
wald( fit, ':Sex') # but no evidence of interaction with ses
wald( fit, '\\\\^2') # nor of curvature

# but we continue for the sake of illustration

L <- Lform( fit, list( 0, 1, 2*ses, 0, Sex == 'Male', (Sex == 'Male')*2*ses), hs)
L
(ww <- wald ( fit, L ))
wald.dd <- as.data.frame(ww, se = 2)
head( wald.dd )

require(lattice)
xyplot( coef + U2 + L2 ~ ses | Sex, wald.dd,
 main= 'Increase in predicted mathach per unit increase in ses')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
