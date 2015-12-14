# gicc
# Author: G. Monette  georges@yorku.ca
# Date: January 3, 2011
# Package: spida



#' Generalized intraclass correlation coefficient for multilevel data %%
#' ~~function to do ... ~~
#' 
#' \code{gicc} facilitates the identification of between-subject variables
#' versus balanced within-subject variables with hierarchical data. Methods
#' handle factors, character variables and numerical variables. For categorical
#' variables, Goodman-Kruskal's tau is returned and for numerical variables,
#' the simple ICC with (variance between) / (variance between + variance
#' within) applied to the rank (by default) of the variable. In either case, a
#' value of 1 signifies a variable that is constant within clusters and a value
#' of 0, a variable that is perfectly balanced within clusters. %% ~~ A concise
#' (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @aliases gicc gicc.default gicc.data.frame gicc.factor
#' @param x a data frame, factor, character variable or a numerical variable %%
#' ~~Describe \code{x} here~~
#' @param by if \code{x} is a data frame, \code{by} can be a formula (e.g.
#' \code{~id}) evaluated in \code{x}. Otherwise, \code{by} is a variable
#' defining clusters. %% ~~Describe \code{by} here~~
#' @param method a character string indicating whether to work with the rank of
#' the raw variable, \code{x}, stripped of missing values, with the raw
#' variable itself, or with \code{is.na(x)}. Can be one of "rank", raw" or
#' "na".
#' 
#' %% ~~Describe \code{method} here~~
#' @param \dots -- not used %% ~~Describe \code{\dots} here~~
#' @return a measure of relative variability within clusters so that 1
#' represents no variability and 0 perfect balance. %% ~Describe the value
#' returned %% If it is a LIST, use %% \item{comp1 }{Description of 'comp1'} %%
#' \item{comp2 }{Description of 'comp2'} %% ...
#' @note %% ~~further notes~~
#' @author G. Monette <georges@@yorku.ca>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#'   gicc( hs, ~ school)
#'   gicc( hs, ~ school, method = 'raw')
#' 
#' @export
gicc <- function(x, by, method,  ...) {

  # Computes a 'generalized intra-class correlation coefficient' using
  # for factors: Goodman and Kruskal's tau based on proportional reduction in
  #    prediction error probability from a proportional prediction rule, and
  # for numerical variables: intraclass correlation coefficient(1,1) based
  #    on ranks of x
  #
  # Args:
  #   x: a variable or a data frame
  #   by: a variable to form cluster, or if x is a data frame, a formula
  #   method: not yet implemented
  # Returns:
  #   The gicc for a variable or for each variable in a data frame
  # TODO: implement different variation indices to see which discriminate
  #       most usefully (e.g. see p. 25 of Agresti (1990) Categorical Data Analysis, Wiley.
    UseMethod("gicc")
}


#' @export
gicc.factor <- function( x, by, method , ...) {
  # See generic function for description
   Verppr <- function( x ) {    # probability of error with proportional prediction rule
      px <- table( x ) / length(x)
      sum( px * (1-px))
   }
   if (method == 'na') { 
     x <- is.na(x)
   } else {
     drop <- is.na(x) | is.na(by)
     x <- x[!drop]
     by <- by[!drop]
   }
   xc <- split(x,by)
   marg.perror <- Verppr(x)
   pcond <- sapply(xc,length)/length(x)
   cond.perror <- sum( pcond * sapply( xc, Verppr))
   #list( marginal = marg.perror, conditional = cond.perror, tau = (marg.perror - cond.perror)/marg.perror)
   (marg.perror - cond.perror)/marg.perror
}

#' @export
gicc.default <- function(x, by, method = 'rank', ...) {
# intraclass correlation based on ranks to reduce influence of outliers
# perhaps we should use normal quantiles
   if ( method == 'na')
      return( gicc( as.factor(x), by, method = 'na' ))
   if (is.character(x))
      return(gicc(as.factor(x), by, method = method))
   drop <- is.na(x) | is.na(by)
   x <-  if(method == "rank") rank(x[!drop]) else x[!drop]
   by <- by[!drop]
   xc <- split(x,by)
   var.bet <- var( sapply( xc, mean))
   wt <- sapply(xc, length) -1
   wt <- wt / sum(wt)
   var.within <- sum(wt * sapply(xc, var), na.rm = TRUE)
   #disp(sapply(xc,var))
   #list(between = var.bet, within = var.within, icc = (var.bet)/(var.bet+var.within))
   var.bet/(var.bet+var.within)
}

#' @export
gicc.data.frame <- function(x, by, method = 'rank', ...) {
# intraclass correlation based on ranks to reduce influence of outliers
# perhaps we should use normal quantiles
  if (inherits(by, "formula"))
    by <- model.frame(by, x, na.action = na.include)
    sapply( x, gicc, by, method)
}
