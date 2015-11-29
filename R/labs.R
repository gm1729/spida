###
###  Row.names heading on data.frame
###

#' Print a labeled data frame
#' 

#' 
#' @param x 
#' @param \dots 
#' @param digits 
#' @param quote 
#' @param right
#' @export 
print.data.frame.lab <-    
  function (x, ..., digits = NULL, quote = FALSE, right = TRUE) 
  {
    labs <- attributes(x)$labs
    if (length(x) == 0) {
      cat("NULL data frame with", length(row.names(x)), "rows\n")
    }
    else if (length(row.names(x)) == 0) {
      print.default(names(x), quote = FALSE)
      cat("<0 rows> (or 0-length row.names)\n")
    }
    else {
      mat <- as.matrix(format.data.frame(x, digits = digits, 
                                         na.encode = FALSE))
      labs <- c(labs,"","")
      labs <- labs[1:2]
      names(dimnames(mat)) <- labs    
      print(mat , ..., quote = quote, right = right)
    }
    invisible(x)
  }

#' @export
"[.data.frame.lab" <- function(x, ...){
  lab <- labs(x)
  ret <- get("[.data.frame")(x,...)
  if( inherits(ret, "data.frame")) labs(ret) <- lab
  ret
} 

#' Add labels -- assignment
#' 
#' @param x 
#' @param \dots 
#' @export
"labs<-" <- function(x,...) UseMethod("labs<-")



#' Add labels -- assignment to data frame
#' 
#' @param x 
#' @param value 
#' @export
"labs<-.data.frame" <- function( x, value ) {
  value <- c( value, "", "") [ 1:2 ]
  attr(x,"labs") <- value
  if( !inherits(x,"data.frame.lab")) class(x) <- c( "data.frame.lab", class(x))
  x
} 
#' Add labels -- default
#' 
#' @param x
#' @param value
#' @export
"labs<-.default" <- function(x, value) {
  nd <- length(dim(x))
  value <- c( value, rep("",nd))[1:nd]
  names(dimnames(x)) <- value
  x
}

#' @export
labs <- function(x,...) UseMethod("labs")


#' Add labels -- data frame
#' 
#' @param x 
#' @param \dots
#' @export
labs.data.frame.lab <- function( x ,...) attr(x,"labs")


#' Add labels -- default
#' 
#' @param x 
#' @param \dots 
#' @export
labs.default <- function(x,...) names(dimnames(x))
