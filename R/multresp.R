##
##  Generating an indicator matrix for 
##    a multiple response variable
##  October 29, 2012
##

##
## create a data frame with indicators according to the presence of a string in a character variable
##

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