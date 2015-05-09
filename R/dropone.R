#
# dropone:
# drop one observation or one cluster at a time
# Value: estimated coeffecients
dropone <- function(fit, ... ) UseMethod('dropone')
dropone.lme <- function (fit, by, verbose = FALSE, ...) 
{
  getSummary <- function( x ){
    re <- getVarCov(x)
    #       disp(re)
    rn <- expand.grid(dimnames(re))
    #       disp(rn)
    rn <- array( apply(rn,1,paste, collapse = "+"), dim = dim(re))
    re <- re[col(re)<=row(re)]
    names( re ) <- rn[col(rn)<=row(rn)]
    ff <- fixef(x)
    names(ff) <- paste('b',names(ff), sep= '.')
    c(ff, dvcov = det(vcov(x)),re, sigma2 = fit$sigma^2 )
  }  
  data <- getData(fit)
  if( missing(by)) by <- 1:nrow(data)
  if (inherits(by, "formula")){ 
    by <-model.frame(by, data, na.action = na.include)
    by <- apply( by, 1, paste, collapse = '/')
  }
  
  by <- as.character(by)
  levs <- unique(by)
  data$.drop <- by
  data.ret <- up( data, ~.drop)
  names(levs) <- levs
  ret <- lapply( levs,
                 function( x) {
                   try(getSummary( update( fit, data = data[by!=x,])))
                 })
  ret <- lapply( ret, function(x) if( inherits(x,"try-error")) NA else x)
  ret <- do.call( rbind,ret)
  disp(dim(ret))
  ret <- as.data.frame(ret)
  ret$.drop <- levs
  disp(dim(ret))
  disp(dim(data.ret))
  ret <- merge(data.ret, ret, by = ".drop")
}


# pairs(z[,grep("^b",names(z))] )