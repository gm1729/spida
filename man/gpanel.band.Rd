\name{gpanel.band}
\alias{gpanel.band}
\alias{gpanel.fit}
\alias{gpanel.labels}
\alias{fillin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Panel functions for predicted values and SE bands
%%  ~~function to do ... ~~
}
\description{
Panel functions for predicted values and SE bands using 'layer' and 'glayer' in the package latticeExtra
}
\usage{
panel.fit(x, y, fit,
           subscripts, ..., type = 'l', font, fontface)

panel.fit(...) # when called within 'layer' or 'glayer' in 'latticeExtra'

panel.band(x, y, lower, upper, subscripts, ..., 
           col, group.number, alpha, 
           col.symbol, border = F, font, fontface) # see examples belor

panel.band(...) # when called within 'layer' or 'glayer' in 'latticeExtra'

panel.labels(x, y, labels , subscripts, ...)

panel.labels(...) # when called within 'layer' or 'glayer' in 'latticeExtra'

fillin(data, form, n = 200, xpd = 1.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
 data frame to be used to add additional values of numeric variable
}
  \item{form}{ formula evaluated in data. The first term defines the variable with values to be filled in and the remaining terms define the variables to be used for grouping determining the minima and maxima within which values are added.
%%     ~~Describe \code{x} here~~
}
\item{n}{ the number of values to be added between the global mininum and maximum. Values falling outside conditional minima and maxima are culled. Default 200.
}
\item{xpd}{
expansion factor to add points beyond minima and maxima. Default 1.0.
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{fit}{ fitted values of a model, generally passed through 'layer' from a call to 'xyplot': e.g. \code{xyplot( y ~ x, data, groups = g, fit = data$yhat, lower = with(data, yhat - 2*se), upper = with(data, yhat + 2*se), subscripts = T)}
%%     ~~Describe \code{lower} here~~
}
  \item{lower}{ 
%%     ~~Describe \code{lower} here~~
}
  \item{upper}{
%%     ~~Describe \code{upper} here~~
}
  \item{subscripts}{ 
%%     ~~Describe \code{subscripts} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{group.number}{
%%     ~~Describe \code{group.number} here~~
}
  \item{alpha}{ to deal with the problem that 'trellis.par.get' will pass a value of alpha = 1, 'panel.band' will set alpha = .3 if alpha = 1. To get a nearly opaque band, use alpha = .99.
%%     ~~Describe \code{alpha} here~~
}
  \item{col.symbol}{ is used to control color when using 'groups'
%%     ~~Describe \code{col.symbol} here~~
}
  \item{border}{ default = FALSE for panel.band.
%%     ~~Describe \code{border} here~~
}
  \item{font}{
%%     ~~Describe \code{font} here~~
}
  \item{fontface}{
%%     ~~Describe \code{fontface} here~~
}
}
\details{
With 'layer' and 'glayer' in 'latticeExtra', these functions can be used to easily generate fitted values and confidence or prediction bands that have a reasonable appearance whether a plot uses 'groups' or not.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The 'panel.bands', 'panel.fit', and 'panel.labels' functions are invoked for their graphical effect.
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

  library(spida)
  library(latticeExtra)
  library(car)
  
  fit <- lm(prestige ~ (income + I(income^2)) * type, Prestige, 
      na.action = na.exclude)
  pred <- cbind( Prestige, predict(fit, newdata = Prestige, se = TRUE))
  head(pred)
  (p <- xyplot( prestige ~ income , pred,  groups = type,
                subscripts = T,
                fit = pred$fit,
                lower = with(pred, fit - 2*se.fit), 
                upper = with(pred, fit + 2*se.fit)))
  p + glayer(gpanel.band(...))
  p + glayer(gpanel.band(...)) + glayer(gpanel.fit(...))
  
  # Another example 
  # to illustrate 'fillin'
  
  fit <- lm( income ~ 
        (education+I(education^2)+I(education^3) +I(education^4))* type, 
        Prestige, na.action = na.exclude)    # overfitting!
  # adding extra values of predictor to get smooth line
  Prestige$occupation <- rownames(Prestige)
  z <- fillin(Prestige, ~education + type,xpd = 1.1)
  
  dim(z)  
  dim(Prestige)
  z <- cbind(z, predict(fit, newdata = z, se = TRUE))
  head(z)
  gd(3,cex=2,lwd=2, alpha = .7)
  (p <-  xyplot( income ~ education, z, groups = type,
                 subscripts = T,
                 fit = z$fit,
                 lower = z$fit - z$se,
                 upper = z$fit + z$se,
                 auto.key = list(space='right', lines= T)))
  p + glayer( gpanel.fit(...))
  p + glayer( gpanel.fit(...)) + glayer(gpanel.band(...))  
  p + glayer( gpanel.band(...,  alpha = .1))
  gd(3,lty=1,lwd=2)
  p + glayer( gpanel.band(...,alpha = .1)) + glayer(gpanel.fit(...))
  
  # with panels and no groups:
  
  (p <-  xyplot( income ~ education| type, z,
                 subscripts = T,
                 fit = z$fit,
                 lower = z$fit - z$se,
                 upper = z$fit + z$se,
                 auto.key = list(space='right', lines= T)))
  p + layer( gpanel.fit(...)) + layer( gpanel.band(...))
  gd(basecol = 'tomato4')
  p + layer( gpanel.band(...,  col = 'grey10')) + layer(gpanel.fit(...))
  p + layer( gpanel.band(...)) + layer(gpanel.fit(...))
  
  # With panels and groups
  
  z <- Prestige 
  z$gender <- with(z, cut( women, c(-1,15,50,101),labels = c("Male","Mixed","Female")))
  tab(z, ~ gender + type)
  z <- fillin( z, ~ education + type + gender, xpd = 1.1)
  fit <- lm( income ~ (education+I(education^2)+I(education^3) )* type * gender, 
             z, na.action = na.exclude)    # overfitting!
  summary(fit)
  z <- cbind( z, predict(fit, newdata = z, se = TRUE))
  head(z)
  (p <-  xyplot( income ~ education| gender, z, groups = type, 
                 subscripts = T,
                 fit = z$fit,
                 lower = z$fit - z$se,
                 upper = z$fit + z$se,
                 layout = c(1,3),
                 auto.key = list(space='right', lines= T, cex = 1.5)))
  
  p + glayer( gpanel.band(...)) + glayer(gpanel.fit(...))
  # trellis.focus()
  # panel.identify(labels= z$occupation)
  # trellis.unfocus()
  
  z$type2 <- with( z, reorder(type,education, mean, na.rm=T))
  gd(3)
  (p <-  xyplot( income ~ education| type2, z, groups = gender, 
                 subscripts = T,
                 fit = z$fit,
                 lower = z$fit - z$se,
                 upper = z$fit + z$se,
                 layout = c(1,3),
                 par.strip.text = list(cex = 2),  
                 auto.key = list(space='right', lines= T, cex = 1.5)))
  
  p + glayer( gpanel.fit(...))
  p + glayer( gpanel.band(...))
  p + glayer( gpanel.band(...)) + glayer(gpanel.fit(...))
  # trellis.focus()
  # panel.identify(labels= z$occupation)
  # trellis.unfocus()
  
  # With panels^2
  # need to remove 'col = col.line'
  
  z <- Prestige 
  z$occ <- rownames(Prestige)
  z$gender <- with(z, cut( women, c(-1,15,50,101),labels = c("Male","Mixed","Female")))
  z$type2 <- with( z, reorder(type,education, mean, na.rm=T))
  tab(z, ~ gender + type2)
  z <- fillin( z, ~ education + type + gender, xpd = 1.1)
  fit <- lm( income ~ (education+I(education^2)+I(education^3) )* type * gender, 
             z, na.action = na.exclude)    # overfitting!
  summary(fit)
  z <- cbind( z, predict(fit, newdata = z, se = TRUE))
  head(z)
  (p <-  xyplot( income ~ education| gender*type, z, 
                 subscripts = T,
                 fit = z$fit,
                 labels = z$occ,
                 lower = z$fit - z$se,
                 upper = z$fit + z$se,
                 auto.key = list(space='right', lines= T, cex = 1.5)))
  
  p + layer( gpanel.fit(...)) + layer( gpanel.band(...))
  p + layer( gpanel.band(..., col = 'black', alpha = .1)) + layer(gpanel.fit(...))  +
    layer(gpanel.text(...))
  

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
