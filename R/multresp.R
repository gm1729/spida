##
##  Generating an indicator matrix for 
##    a multiple response variable
##  October 29, 2012
##

##
## create a data frame with indicators according to the presence of a string in a character variable
##
#' Create an indicator matrix from a multiple response variable
#' 
#' A multiple reponse variable is one whose values consist of a list of
#' possible choices. Such a variable can be represented in a linear model
#' through a matrix of indicator variables.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param x %% ~~Describe \code{x} here~~ a factor or character vector
#' containing comma- or blank-separated lists of reponse choices for each
#' subject. For example, each subject could have a list of symptoms: c('cough
#' fever', 'cough bronchitis', 'headache wheeze', '', 'cough high_fever', NA).
#' @param split %% ~~Describe \code{split} here~~ a regular expression matching
#' the separator between fields. The default is \code{"[ \t,]+"} which matches
#' one more blank, tab or comma.  To use semicolons instead of commas as
#' separators, use: \code{"[ \t;]+"} .
#' @param sep the separator used to generated the column names of the indicator
#' matrix.  \code{sep} separates the root from the response values. Default:
#' \code{"."}. For example, if \code{x} is the name of the variable, then the
#' default column name for response value 'cough' is \code{x.cough}.
#' @param root the root for column names. By default, the name of the variable.
#' the
#' @return A matrix of indicator variables. For, example, if \code{x <-c('cough
#' fever', 'cough bronchitis', 'headache wheeze', '', 'cough high_fever', NA)}
#' \code{multresp(x)} produces: \code{ x.bronchitis x.cough x.fever x.headache
#' x.high_fever x.wheeze [1,] 0 1 1 0 0 0 [2,] 1 1 0 0 0 0 [3,] 0 0 0 1 0 1
#' [4,] 0 0 0 0 0 0 [5,] 0 1 0 0 1 0 [6,] 0 0 0 0 0 0
#' } 
#' @examples
#' 
#' x <-c('cough fever', 'cough bronchitis', 'headache wheeze', '', 'cough high_fever', NA)
#' multresp(x)
#' 
#' data <- data.frame(x = x, y = 1:6)
#' data <- with(data, cbind(data, multresp(x)))
#' data
#' 
#' @export
multresp <- function(x,split="[ \t,]+",sep='.',root=deparse(substitute(x))) {
  root
  x <- as.character(x)
  sp <- strsplit(x,split=split)
  responses <- na.omit(unique(unlist(sp)))
  nams <- make.names(paste(root,sep,responses,sep=""), unique = TRUE)
  suppressWarnings( responses.num <- as.numeric(responses))
  n.num <- sum(!is.na(responses.num))
  n.char <- sum(is.na(responses.num))
  responses.char <- responses
  responses.char[!is.na(responses.num)] <- NA
  ord.num <- order(responses.num)
  ord.char <- order(responses.char)
  ord <- c(ord.num[seq_len(n.num)],ord.char[seq_len(n.char)])
  xactmatch <- function(pat,patlist) {
    # returns TRUE if exact match in any element of patlist
    sapply(patlist, function (x) any(sapply( x, identical, pat)))
  }
  mat <- 1*sapply( responses, xactmatch, sp)
  colnames(mat) <- nams
  mat[,ord]
}
