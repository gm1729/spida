###
### Copyright 2012  Georges Monette <georges@yorku.ca>
###
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

###
### Lfx.R
### part of the spida package
#
# 2012 12 05:
#  Modelled and improved from Lform: Allows estimation of effects of arbitrary
#  order evaluated over a data frame.  Correctly handles interaction terms involving
#  matrices. 
#  Lfx generates a matrix that can be used as a linear hypothesis matrix in 'wald'
getData <- function(x,...) UseMethod("getData")

Lfx <-     
  function (fit, expr.list, data = getData(fit), prefix = "1 * ", wrap = FALSE, debug = FALSE) 
  {
    # 2012-12-03:
    # 1. modified to work with M with factor instead of formula
    # 2. Lform(fit) includes comments from Lform0
    # 3. More informative error messages 
    
    Lfx0 <- function(fit, factors = getFactorNames(fit), wrap = FALSE,
                       debug = debug) {
      ts <- colnames(attr(terms(fit),'factors'))
      ts <- strsplit(ts,":")
      if(debug) disp(ts)
      # add functions as well because some might generate matrices:
      
      facl <- lapply( ts, function(x) x %in% factors | grepl(")",x))
      if(debug) disp(facl)
      ret <- lapply( seq_along(facl),  function(i){
        lapply(seq_along(facl[[i]]),
               function(j) if(facl[[i]][j]) 
                 paste("M(",ts[[i]][j],")", sep = "") else ts[[i]][j] 
        )})
      
      ret <- lapply( ret, paste, collapse = " * ")    
      ret <- lapply( ret, function(s) paste(prefix,s,sep=""))  # add prefix 
      ret <- do.call( paste, c("list( 1",ret, sep = ",\n"))
      ret <- paste(ret,"\n)")
      if (wrap) ret <- paste("with(data, do.call( cbind,\n",
                             ret, "\n ))\n")
      cat(ret)
      invisible(ret)
    }
    
    
    if (missing(expr.list)) return( Lfx0(fit, wrap = wrap, debug = debug))
    gg <- getFix(fit)
    Lsub <- do.call(cbind, eval(substitute(expr.list), data))
    L <- Lsub
    rownames(L) <- rownames(data)
    colnames(L) <- names(gg$fixed)
    
    attr(L, "data") <- data
    L
  }

# 2012 12 03: Replaced with following functions
# M <- function(form, ..., keep.intercept = FALSE) {
#   ret <- model.matrix(form,...)
#   if( !keep.intercept ){
#     pos <- grep( '^\\(Intercept\\)$', colnames(ret))
#     ret <- if( length(pos) ) ret[, -pos] 
#   }
#   ret
# }

# To do: return object of type "M" and
# define '*' method to do column by column multiplication
# -- DONE

M <- function(x,...) UseMethod("M")

M.formula <- function (form, ..., keep.intercept = FALSE) { 
  ret <- model.matrix(form, ...)
  if (!keep.intercept) {
    pos <- grep("^\\(Intercept\\)$", colnames(ret))
    ret <- if (length(pos)) 
      ret[, -pos]
  }
  class(ret) <- 'M'
  ret
}

M.factor <- function( x, base = NULL) {
  ret <- contrasts(x)[x,]
  if( ! is.null(base)) {
    if ( !is.factor(base)) stop("base must be a factor")
    if (!identical(levels(x), levels(base))) warning("'x' and 'base' should have identical levels")
    ret <- ret - contrasts(base)[base,]
  }
  class(ret) <- 'M'
  ret
}

M.M <- function(x, ...) x

M.default <- function(x,...) {
  x <- cbind(x)
  class(x) <- c('M',attr(x,"class"))
  x
}

"*.M" <- function (e1, e2) {
  # column by column multiplication
  # to expand interaction terms
  # LHS runs fastest
  if( !inherits(e1,"M") && is.numeric(e1) && length(e1) == 1) {
    return( M( e1 * unclass(e2)))    
  }
  if( !inherits(e2,"M") && is.numeric(e2) && length(e2) == 1) {
    return( M( unclass(e1) * e2))    
  }
  e1 <- cbind(e1)
  e2 <- cbind(e2)
  n1 <- ncol(e1)
  n2 <- ncol(e2)
  ret <- c(e1) * c(apply(e2, 2, rep, n1))
  ret <- matrix(ret, nrow = nrow(e1))
  #  rownames(ret) <- rownames(e1)
  if (!is.null(c1 <- colnames(e1)) && !is.null(c2 <- colnames(e2))) {
    cn <- paste(rep(c1, length(c2)), rep(c2, each = length(c1)), 
                sep = ":")
    colnames(ret) <- cn
  }
  class(ret) <- "M"
  ret
}

# To make is easier to generate pairwise comparisons without
# redundancies
# e.g.
#  > pred <- expand.grid( type = levels(Prestige$type),
#                        type0 = levels(Prestige$type),
#                        education = seq(6,18,3))
#  > pred <- subset( pred, type0 < type)


"<.factor" <- function(x,y) {
  if ( ! identical( levels(x), levels(y))) warning( "Different 'levels' in x and y")
  as.numeric(x) < as.numeric(y)
}


">.factor" <- function(x,y) {
  if ( ! identical( levels(x), levels(y))) warning( "Different 'levels' in x and y")
  as.numeric(x) > as.numeric(y)
}

"<=.factor" <- function(x,y) {
  if ( ! identical( levels(x), levels(y))) warning( "Different 'levels' in x and y")
  as.numeric(x) <= as.numeric(y)
}

">=.factor" <- function(x,y) {
  if ( ! identical( levels(x), levels(y))) warning( "Different 'levels' in x and y")
  as.numeric(x) >= as.numeric(y)
}

