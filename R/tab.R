###
###
### From spida-test
### Modified by GM 2013-08-20
### 1. Added ability to use pr and pct in Tab
###    dropping total margins but keeping by default "All"
###    margins
### 2. Added utility function 'dropLast'
### Modified by GM 2009-09-26
### tab( ..., pct = 1) standardizes rows so thay add to 100
### tab( ..., pr = 1) standardizes rows so thay add to 1
### tab( ..., pr = 0) standardizes entire table so it sums to 1
###
###    Incorporate testing from package rdc
####

# Changes:
# 2014:
#   October 16: added tab.table, tab.matrix, tab.array
#
# TODO: update tab.Rd to explain new features in Tab
# and keep = "All"



#' @export
.mat2arr <- function(x) {
      ret <- as.list(x)
      dim(ret) <- NULL
      nams <- expand.grid( dimnames(x))
      for ( ii in 1:length(nams)) {
          nams[[ii]] <- paste( names(nams)[ii] , nams[[ii]], sep = " = ")
      }
      nams <- c( nams, sep = ", ")
      nams <- do.call( paste, nams)
      names(ret) <- nams
      ret
}
#mat2arr( zza)


#' @export
dropLast <- function( mat ,drop = FALSE, keep = NULL) {
  # NEW: 2013-08-20, G. Monette
  # utility function to drop totals from tab
  ind.last <- dim(mat)
  ind.keep <- as.list(- ind.last) # drop last index of each dim unless name in keep
  last <- function(x) x[length(x)]
  names.last <- sapply( dimnames(mat), last)
  if ( !is.null(keep)) {
    no.drop <- names.last %in% keep
    ind.keep[no.drop] <- lapply( ind.last[no.drop], seq_len)  
  }
  call <- c(list(mat),ind.keep,drop=drop)
  # print(call)
  do.call( `[`,call)
}

#' @export
Tab <- function(...,keep = "All") {
  # New version of Tab that handles pct and pr
  # To keep the "All" and not the "Total" rows,
  # specify keep = "All"
  # BUGS: would be more efficient if it 
  #       called tab(...,total.margins=FALSE)
  #       when pct or pr arguments are not given
  dropLast(tab(...), keep = keep)
}

#' @export
pab <- Tab     # legacy




#' Table of frequencies or relative frequencies bordered with totals and
#' including NAs
#' 
#' Generates a table of frequencies or relative frequencies or relative
#' percentages %% ~~ A concise (1-5 lines) description of what the function
#' does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @aliases tab tab.formula tab.data.frame tab.default Tab
#' @param \dots as with \code{table}, one or more objects which can be
#' interpreted as factors (including character strings), or a list (or data
#' frame) whose components can be so interpreted.
#' @param data a data frame in which formula are interpreted
#' @param fmla a formula whose right-hand side names the variables to be used
#' for tabulation. The optional left-hand side specifies a variable to be used
#' for weights.
#' @param useNA whether to include NA levels. The default is "ifany". Can also
#' be set to "no" or "always".
#' @param pct margins to be scaled to sum to 100
#' @param pr margins to be scaled to sum to 1
#' @param total.margins if FALSE, generate table without margins
#' @param weights instead of generating a frequency table, generate a table
#' with the sum of the weights
#' @param keep names of margins to keep with 'Tab', default = "All". To drop
#' all margins, use default = "".
#' @return An object of class 'table' of dimension equal to the number of
#' variables, with optional margins showing totals. Elements of the matrix can
#' be frequencies, relative frequencies, percentages or sums of weights.
#' @note %% ~~further notes~~
#' @author Georges Monette
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' titanic <- as.data.frame(Titanic)
#' head(titanic)
#' tab(titanic, ~ Sex + Age, weights = Freq) 
#' 
#' tab(titanic, ~ Sex + Survived + Age, weights = Freq)
#' 
#' round(tab(titanic, ~ Sex + Survived + Age, 
#'     pct = c(1,3), weights = Freq),2)
#' 
#' round(Tab(titanic, ~ Sex + Survived + Age, 
#'     pct = c(1,3), weights = Freq),2)
#' 
#' round(Tab(titanic, ~ Sex + Survived + Age, 
#'     pct = c(1,3), weights = Freq, keep = ""),2)
#' 
#' @export
tab <- function(x,...) UseMethod("tab")

#' @export
tab.table <- function(x,...) {
    data <- as.data.frame(x)
    wt <- data$Freq
    data$Freq <- NULL
    tab(data, ..., weights = wt)
}

#' @export
tab.matrix <- function(x,...) tab(as.table(x),...)
#' @export
tab.array <- function(x,...) tab(as.table(x),...)


#' @export
tab.formula <- function( fmla, data = sys.frame(sys.parent()), ... ) tab.data.frame(data,fmla,...)

#' @export
tab.data.frame <- function (dd, fmla, 
                            total.margins = TRUE, 
                            useNA = "ifany",  
                            pct = NULL, pr = NULL, 
                            test = FALSE,
                            weights = NULL,
                            na.rm = NULL,
                            all.label = "All")
{
  # GM: 2014 08 22: modified handling of lhs to fix bug when variable in lhs
  #                 also appears in rhs
  if (missing(fmla)) {
    fmla <- parse(text = paste(c("~",paste(names(dd),collapse="+"))))
    fmla <- eval(fmla)
    environment(fmla) <- parent.frame()
  }
  if (is.null(weights) && (length(fmla) >2 )) {
    weights <- model.frame(fmla[-3], dd, na.action = NULL)[[1]]
    xx <- model.frame(fmla[-2], dd, na.action = NULL)
  } else {
    xx = model.frame(fmla, dd, na.action = NULL)
    weights <- eval(substitute(weights), dd, environment(fmla))    
  }
  if(!is.null(weights) && any(is.na(weights))) {
    warning("NAs in weights set to -Inf")
    weights[ is.na(weights) ] <- -Inf
  }
  #  disp(weights)
  xx = c(xx, list( total.margins = total.margins, useNA = useNA, 
                   pct=pct, pr = pr, test=test, na.rm = na.rm,
                   weights = weights,all.label=all.label) )
  do.call("tab", xx)
}

#' @export
tab.default <- function (..., total.margins = TRUE, 
                         pct = NULL, pr = NULL, 
                         useNA = "ifany",
                         test = FALSE, 
                         weights = NULL,
                         na.rm = NULL,
                         all.label = "All")
{
  
  if(!is.null(na.rm)) useNA <- if(na.rm) "no" else "ifany" 
  aa <- list(...)
  if (length(aa) == 1 && is.list(aa[[1]])) {
    aa <- c(aa[[1]],list( total.margins = total.margins,  
                          useNA = useNA, pr = pr, pct = pct, test=test,
                          na.rm = na.rm, weights = weights))
    disp(aa[[1]])
    return(do.call("tab", aa[[1]]))
  }
  if (is.null(names(aa))) {
    nns = names(match.call())
    names(aa) = nns[2:(1 + length(aa))]
  }
  if(useNA=="ifany") for (ii in 1:length(aa)) aa[[ii]] <- factor(aa[[ii]], exclude = NULL)  # table uses 'ifany' correctly when a number but not for factors
  if( is.null(weights)){
    aa[["useNA"]] <- useNA 
    aa[["na.rm"]] <- NULL
    ret <- do.call("table", aa)
  } else {
    ret <- as.table(tapply(weights, aa, sum ) )
    ret[is.na(ret)] <- 0
  }
  if (test) {
    if( length(dim(ret)) < 3) {
      test.out <- chisq.test(ret)
    } else if ( length(dim(ret)) == 3)  {
      test.out <- apply( ret, 3:length(dim(ret)), chisq.test)
      
    } else if ( length(dim(ret)) > 3)  {
      test.out <- .mat2arr(apply( ret, 3:length(dim(ret)), chisq.test))
      
    }
  }
  if ( !is.null(pr)) ret <- acond( ret, MARGIN = pr, all.label = all.label)
  else if ( !is.null(pct)) ret <- 100* acond( ret, MARGIN = pct, all.label = all.label)
  else if (total.margins) ret = atotal(ret)
  if( test ) attr(ret,'test') <- test.out
  if( test ) ret <- unclass(ret)   # so attributes will print  
  ret
}





#' Transform a frequency table into relative frequencies relative to a margin.
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param x %% ~~Describe \code{x} here~~
#' @param MARGIN %% ~~Describe \code{MARGIN} here~~
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' ## The function is currently defined as
#' function (x, MARGIN = NULL) 
#' {
#'     debug <- F
#'     x <- atotal(x)
#'     if (is.null(MARGIN)) 
#'         return(x/x[length(x)])
#'     d <- dim(x)
#'     if (debug) {
#'         cat("dim(x)\n")
#'         print(d)
#'     }
#'     dn <- dimnames(x)
#'     if (debug) {
#'         cat("dimnames(x)\n")
#'         print(dn)
#'     }
#'     n <- length(d)
#'     m <- length(MARGIN)
#'     perm <- c(MARGIN, setdiff(1:n, MARGIN))
#'     perm.inverse <- order(perm)
#'     if (debug) 
#'         print(perm)
#'     if (debug) 
#'         print(perm.inverse)
#'     x <- aperm(x, perm)
#'     d <- dim(x)
#'     zl <- list(x)
#'     for (ii in 1:length(d)) zl[[ii + 1]] <- seq(if (ii <= m) 
#'         1
#'     else d[ii], d[ii])
#'     tots <- do.call("[", zl)
#'     ret <- array(c(x)/c(tots), dim = d)
#'     ret <- aperm(ret, perm.inverse)
#'     if (debug) 
#'         print(dim(ret))
#'     if (debug) 
#'         print(dn)
#'     dimnames(ret) <- dn
#'     ret
#'   }
#' 
#' @export
acond <- function (x, MARGIN = NULL, total.margins = TRUE, all.label = "All")
{
    debug <- F
    x <- if( total.margins ) atotal(x)
    if (is.null(MARGIN)|| max(MARGIN) < 1)
        return(x/x[length(x)])
    d <- dim(x)
    if (debug) {
        cat("dim(x)\n")
        print(d)
    }
    dn <- dimnames(x)
    if (debug) {
        cat("dimnames(x)\n")
        print(dn)
    }
    n <- length(d)
    m <- length(MARGIN)
    if ( length(setdiff( MARGIN, 1:n)) > 0) stop("MARGIN must select dimensions of x, or be equal to 0")
    perm <- c(MARGIN, setdiff(1:n, MARGIN))
    perm.inverse <- order(perm)
    if (debug)
        disp(perm)
    if (debug)
        disp(perm.inverse)
    x <- aperm(x, perm)
    d <- dim(x)
    zl <- list(x)
    for (ii in 1:length(d)) zl[[ii + 1]] <- seq(if (ii <= m)
        1
    else d[ii], d[ii])
    tots <- do.call("[", zl)
    ret <- array(c(x)/c(tots), dim = d)
    ret <- aperm(ret, perm.inverse)
    if (debug)
        disp(dim(ret))
    if (debug)
        disp(dn)
    for ( ii in MARGIN) dn[[ii]][length(dn[[ii]])] <- all.label
    dimnames(ret) <- dn
    ret
}

# environment(acond) <- as.environment("package:spida.beta")




#' Relative proportions
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param x %% ~~Describe \code{x} here~~
#' @param MARGIN %% ~~Describe \code{MARGIN} here~~
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' ## The function is currently defined as
#' function (x, MARGIN = NULL) 
#' {
#'     debug <- F
#'     x <- atotal(x)
#'     if (is.null(MARGIN)) 
#'         return(x/x[length(x)])
#'     d <- dim(x)
#'     if (debug) {
#'         cat("dim(x)\n")
#'         print(d)
#'     }
#'     dn <- dimnames(x)
#'     if (debug) {
#'         cat("dimnames(x)\n")
#'         print(dn)
#'     }
#'     n <- length(d)
#'     m <- length(MARGIN)
#'     perm <- c(MARGIN, setdiff(1:n, MARGIN))
#'     perm.inverse <- order(perm)
#'     if (debug) 
#'         print(perm)
#'     if (debug) 
#'         print(perm.inverse)
#'     x <- aperm(x, perm)
#'     d <- dim(x)
#'     zl <- list(x)
#'     for (ii in 1:length(d)) zl[[ii + 1]] <- seq(if (ii <= m) 
#'         1
#'     else d[ii], d[ii])
#'     tots <- do.call("[", zl)
#'     ret <- array(c(x)/c(tots), dim = d)
#'     ret <- aperm(ret, perm.inverse)
#'     if (debug) 
#'         print(dim(ret))
#'     if (debug) 
#'         print(dn)
#'     dimnames(ret) <- dn
#'     ret
#'   }
#' 
#' @export
aprop <- acond    # older name


