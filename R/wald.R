
## wald.R
## This a collection of functions designed to facilitate testing hypotheses
## with Wald tests.
## The methods are readily extended to any fitting method that returns a vector
## of estimated parameters and an estimated variance matrix.
## Extensions are implemented through the 'getFix' generic function.


## see also Leff
## Changes:
## 2014 06 04: changed fit@fixef to fixef(fit) in a number of 'getFix' methods
## October 2, 2011:  modified 'wald' to better handle rank-deficient models
##                   previously columns of L and columns and rows of vcov
##                   corresponding to NAs in beta were deleted ignoring the
##                   possibility that L beta is not estimable if any 
##                   non-zero element of L multiplies an NA in beta.
##
## 2013 06 18: added df argument to wald to override denominator dfs. Useful
##             for saturated binomial fits or other binomials where residual
##             dfs are not appropriate.
##
## 2013 09 17: added getFix.multinom with df = Inf
##                  

glh <- function( ...) wald( ...)    # previous name for 'wald' function

wald <- function(fit, Llist = "",clevel=0.95, data = NULL, debug = FALSE , maxrows = 25,
     full = FALSE, fixed = FALSE, invert = FALSE, method = 'svd',df = NULL) {
        if (full ) return( wald ( fit, model.matrix(fit)))
          dataf <- function(x,...) {
           x <- cbind(x)
           rn <- rownames(x)
           if( length( unique(rn)) < length(rn)) rownames(x) <- NULL
           data.frame(x,...)
         }
     as.dataf <- function(x,...) {
           x <- cbind(x)
           rn <- rownames(x)
           if( length( unique(rn)) < length(rn)) rownames(x) <- NULL
           as.data.frame(x,...)
     }

     unique.rownames <- function(x) {
                ret <- c(tapply(1:length(x), x , function(xx) {
                    if ( length(xx) == 1) ""
                    else 1:length(xx)
                })) [tapply(1:length(x),x)]
                ret <- paste(x,ret,sep="")
                ret
     }
#      if(debug) disp( Llist)
     if( is.character(Llist) ) Llist <- structure(list(Llist),names=Llist)
     if(!is.list(Llist)) Llist <- list(Llist)

     ret <- list()
     fix <- getFix(fit)
#      if(debug) disp(fix)
     beta <- fix$fixed
     vc <- fix$vcov

     dfs <- if(is.null(df) ) fix$df else df + 0*fix$df
#      if(debug) disp(Llist)
     for (ii in 1:length(Llist)) {
         ret[[ii]] <- list()
         Larg <- Llist[[ii]]
#          if(debug) {
#               disp(ii)
#               disp(Larg)
#          }
         L <- NULL
         if ( is.character(Larg)) {
            L <- Lmat(fit,Larg, fixed = fixed, invert = invert)
         } else {
            if ( is.numeric(Larg)) {
               if ( is.null(dim(Larg))) {
                  if(debug) disp(dim(Larg))
                  if ( (length(Larg) < length(beta)) && (all(Larg>0)||all(Larg<0)) ) {
                        L <- diag(length(beta))[Larg,]
                        dimnames(L) <- list( names(beta)[Larg], names(beta))
                  } else L <- rbind( Larg )
               }
               else L <- Larg
            }
         }
#          if (debug) {
#             disp(Larg)
#             disp(L)
#          }
         ## Delete coefficients that are NA
         Ldata <- attr( L , 'data')
         
         ## identify rows of L that are not estimable because they depend on betas that are NA
         Lna <- L[, is.na(beta), drop = FALSE]
         narows <- apply(Lna,1, function(x) sum(abs(x))) > 0
         
         L <- L[, !is.na(beta),drop = FALSE]
         attr(L,'data') <- Ldata
         beta <- beta[ !is.na(beta) ]
         
         ## Anova
         if( method == 'qr' ) {
             qqr <- qr(t(na.omit(L)))     
             #Qqr <- Q(t(L))
             L.rank <- qqr$rank
             #L.rank <- attr(Qqr,'rank')
             #L.miss <- attr(Qqr,'miss')
             if(debug)disp( t( qr.Q(qqr)))
             L.full <- t(qr.Q(qqr))[ 1:L.rank,,drop=FALSE]
             #L.full <- t(Qqr[!L.miss,])[ 1:L.rank,,drop=F]
         } else if ( method == 'svd' ) {
              if(debug) disp(L)
#              if(debug)disp( t(na.omit(t(L))))
#              sv <- svd( t(na.omit(t(L))) , nu = 0 )
             sv <- svd( na.omit(L) , nu = 0 )
 
              if(debug)disp( sv )
             tol.fac <- max( dim(L) ) * max( sv$d )
             if(debug)disp( tol.fac )
             if ( tol.fac > 1e6 ) warning( "Poorly conditioned L matrix, calculated numDF may be incorrect")
             tol <- tol.fac * .Machine$double.eps
             if(debug)disp( tol )
             L.rank <- sum( sv$d > tol )
             if(debug)disp( L.rank )
             if(debug)disp( t(sv$v))
             L.full <- t(sv$v)[seq_len(L.rank),,drop = FALSE]
         } else stop("method not implemented: choose 'svd' or 'qr'")
         
         # from package(corpcor)
         # Note that the definition tol= max(dim(m))*max(D)*.Machine$double.eps
         # is exactly compatible with the conventions used in "Octave" or "Matlab".


         if (debug && method == "qr") {
            disp(qqr)
            disp(dim(L.full))
            disp(dim(vc))
            disp(vc)
         }

         if (debug) disp(L.full)
         if (debug) disp( vc )
         
         vv <-  L.full %*% vc %*% t(L.full)
         eta.hat <- L.full %*% beta
         Fstat <- (t(eta.hat) %*% qr.solve(vv,eta.hat,tol=1e-10)) / L.rank
         included.effects <- apply(L,2,function(x) sum(abs(x),na.rm=TRUE)) != 0
         denDF <- min( dfs[included.effects])
         numDF <- L.rank
         ret[[ii]]$anova <- list(numDF = numDF, denDF = denDF,
                         "F-value" = Fstat,
                         "p-value" = pf(Fstat, numDF, denDF, lower.tail = FALSE))
         ## Estimate
         etahat <- L %*% beta
        # NAs if not estimable:
        etahat[narows] <- NA
         if( nrow(L) <= maxrows ) {
              etavar <- L %*% vc %*% t(L)
              etasd <- sqrt( diag( etavar ))
         } else {
              etavar <- NULL
              etasd <- sqrt( apply( L * (L%*%vc), 1, sum))
         }

         denDF <- apply( L , 1 , function(x,dfs) min( dfs[x!=0]), dfs = dfs)
         
         aod <- cbind( Estimate=c(etahat),
             Std.Error = etasd,
             DF = denDF,
             "t-value" = c(etahat/etasd),
             "p-value" = 2*pt(abs(etahat/etasd), denDF, lower.tail =FALSE))
          colnames(aod)[ncol(aod)] <- 'p-value'
       if (debug ) disp(aod)
       if ( !is.null(clevel) ) {
          #print(aod)
          #print(aod[,'DF'])
          #print(aod[,'etasd'])
           hw <- qt(1 - (1-clevel)/2, aod[,'DF']) * aod[,'Std.Error']
           #print(hw)
           aod <- cbind( aod, LL = aod[,"Estimate"] - hw, UL = aod[,"Estimate"] + hw)
           #print(aod)
           if (debug ) disp(colnames(aod))
           labs <- paste(c("Lower","Upper"), format(clevel))
           colnames(aod) [ ncol(aod) + c(-1,0)] <- labs
       }
       if (debug ) disp(rownames(aod))
       aod <- as.dataf(aod)

       rownames(aod) <- rownames(as.dataf(L))
       labs(aod) <- names(dimnames(L))[1]
       ret[[ii]]$estimate <- aod
       ret[[ii]]$coef <- c(etahat)
       ret[[ii]]$vcov <- etavar
       ret[[ii]]$L <- L
       ret[[ii]]$se <- etasd
       ret[[ii]]$L.full <- L.full
       ret[[ii]]$L.rank <- L.rank
       if( debug ) disp(attr(Larg,'data'))
       ret[[ii]]$data <- attr(Larg,'data')
  }
  names(ret) <- names(Llist)
  attr(ret,"class") <- "wald"
  ret
}

print.wald <- function(x,round = 6, pround = 5,...) {
    pformat <- function(x, digits = pround) {
        x <- format(xx <- round(x,digits))
        x[ as.double(xx) == 0 ] <- paste(c("<.",rep('0',digits-1),'1'),collapse="")
        x
    }
    rnd <- function(x,digits) {
        if (is.numeric(x)) x <- round(x,digits=digits)
        format(x)
    }
             for( ii in 1:length(x)) {
                  nn <- names(x)[ii]
                  tt <- x[[ii]]
                  ta <- tt$anova

                  ta[["p-value"]] <- pformat(ta[["p-value"]])
                  print(as.data.frame(ta,row.names=nn))
                  te <- tt$estimate
                  rowlab <- attr(te,"labs")

                  te[,'p-value'] <- pformat( te[,'p-value'])
                  if ( !is.null(round)) {
                     for ( ii in 1:length(te)) {
                         te[[ii]] <- rnd(te[[ii]],digits=round)
                     }
                  }
                  labs(te) <- rowlab
                  print(te,digits=round,...)
                  cat("\n")
             }
        invisible(x)
}

as.data.frame.wald <- function( obj, se , digits = 3 , sep = "" , which = 1 ) {
# modified by GM 2010_09_20 to avoid problems with coefs with duplicate rownames
         dataf <- function(x,...) {
           x <- cbind(x)
           rn <- rownames(x)
           if( length( unique(rn)) < length(rn)) rownames(x) <- NULL
           data.frame(x,...)
         }
         obj = obj [which]
         if ( length(obj) == 1) {

            cf <- obj[[1]]$coef
            ret <- if ( missing(se)) data.frame( coef = cf, se = obj[[1]]$se)
            else {
               if ( is.null( names(se))) names(se) <-
                             sapply(se,function(x) as.character(round(x,digits)))
               SE <- obj[[1]]$se
               SEmat <- cbind(SE) %*% rbind(se)
               cplus <- cf + SEmat
               cminus <- cf - SEmat
               colnames(cplus) <- paste( "U",colnames(cplus),sep=sep)
               colnames(cminus) <- paste( "L",colnames(cminus),sep=sep)
               ret <- data.frame( coef = cf)
               ret <- cbind( ret, cplus, cminus)
            }
            if( is.null(dd <- obj[[1]]$data)) return( ret)
            else return( cbind(ret, dd))
         }
         else ret <- lapply( obj, as.data.frame.wald)
         ret
}

coef.wald <- function( obj , se = FALSE ) {
         if ( length(obj) == 1) {
            ret <-
            ret <- obj[[1]]$coef
            if ( is.logical(se) && (se == TRUE) ) {
               ret <- cbind( coef = ret, se = obj[[1]]$se)

            } else if ( se > 0 ){
               ret <- cbind( coef = ret, coefp = ret+se*obj[[1]]$se,
                      coefm = ret - se*obj[[1]]$se)
               attr(ret,'factor') <- se
            }
         }
         else ret <- sapply( obj, coef.wald )
         ret
}

##
##
##   Functions to perform a GLH on lm, lme or lmer models
##   August 13, 2005
##
##
##
##   Lmat: generate a hypothesis matrix based on a pattern
##
##   glh
##   Lmat
##   Ldiff
##   getFix
##
##   print.glh
##


getFix <- function(fit,...) UseMethod("getFix")

getFix.multinom <- function(fit,...) {
  ret <- list()
  ret$fixed <- c(t(coef(fit)))
  ret$vcov <- vcov(fit)
  names(ret$fixed) <- rownames(ret$vcov)
  df <- nrow(na.omit(cbind(fit$residuals))) - length(ret$fixed)
  ret$df <- rep( df, length(ret$fixed))
  ret
}


getFix.lm <- function(fit,...) {
       ss <- summary(fit)
       ret <- list()
       ret$fixed <- coef(fit)
       ret$vcov <- ss$sigma^2 * ss$cov.unscaled
       ret$df <- rep(ss$df[2], length(ret$fixed))
       ret
}


getFix.glm <- function(fit,...) {
       ss <- summary(fit)
       ret <- list()
       ret$fixed <- coef(fit)
       ret$vcov <- vcov(fit)
       ret$df <- rep(ss$df.residual, length(ret$fixed))
       ret
}

getFix.lme <- function(fit,...) {
       require(nlme)
       ret <- list()
       ret$fixed <- nlme::fixef(fit)
       ret$vcov <- fit$varFix
       ret$df <- fit$fixDF$X
       ret
}


getFix.gls <- function(fit,...) {
  require(nlme)
  ret <- list()
  ret$fixed <-coef(fit)
  ret$vcov <- vcov(fit)
  ds <- fit$dims
  df <- ds[[1]] - sum( unlist( ds[-1]))
  ret$df <- rep(df, length(coef(fit)))
  ret
}


getFix.lmer <- function(fit,...) {
  # 2014 06 04: changed fit@fixef to fixef(fit)     
  ret <- list()
       ret$fixed <- fixef(fit)
       ret$vcov <- as.matrix( vcov(fit) )
       # ret$df <- Matrix:::getFixDF(fit)
       ret$df <- rep( Inf, length(ret$fixed))
       ret
}

getFix.glmer <- function(fit,...) {
  # 2014 06 04: changed fit@fixef to fixef(fit)     
  
  ret <- list()
       ret$fixed <- fixef(fit)
       ret$vcov <- as.matrix(vcov(fit))
       # ret$df <- Matrix:::getFixDF(fit)
       ret$df <- rep( Inf, length(ret$fixed))
       ret
}


getFix.mer <- function(fit,...) {
  # 2014 06 04: changed fit@fixef to fixef(fit)     
  
       ret <- list()
       ret$fixed <- fixef(fit)
       ret$vcov <- as.matrix(vcov(fit))
       # ret$df <- Matrix:::getFixDF(fit)
       ret$df <- rep( Inf, length(ret$fixed))
       ret
}

getFix.zeroinfl <- function(fit,...){
       ret <- list()
       ret$fixed <- coef(fit)
       ret$vcov <- as.matrix(vcov(fit))
       # ret$df <- Matrix:::getFixDF(fit)
       ret$df <- rep( Inf, length(ret$fixed))
       ret
}

getFix.mipo <- function( fit, ...){ 
  # pooled multiple imputation object in mice
  # uses the minimal df for components with non-zero weights
  # -- this is probably too conservative and should 
  # improved
  ret <- list()
  ret$fixed <- fit$qbar
  ret$vcov <- fit$t
  ret$df <- fit$df
  ret
}






Vcov <- function(fit) {
     getFix(fit)$vcov
}

Vcor <- function(fit) {
     vc <- cov2cor(getFix(fit)$vcov)
     svds <- svd(vc)$d
     attribute(vc,'conditionNumber') <- svds[1]/svds[length(svds)]
     vc
}


###  getFix: function designed to be used internally to get coef, var(coef) and df.resid

getFix.rdc <- function(fit, ...) UseMethod("getFix")

getFix.rdc.lmer <- function(fit, ...) {
            ret <- list()
            ret$fixed <- fixef(fit)
            ret$vcov <- vcov(fit)
            ret$df <- Matrix:::getFixDF(fit)
            ret
}

getFix.rdc.lm <- function(fit, ...) {
          ret <- list()
          ret$fixed <- coef(fit)
          ret$vcov <- vcov(fit)
          ret$df <- fit$df.residuals
          ret
}


getFix.MCMCglmm <- function(fit,...) {
    ret <- list()
    ret$fixed <- apply(fit$Sol, 2, mean)
    ret$vcov <- var( fit $ Sol)
    ret$df <- rep(Inf, length(ret$fixed))
    ret
}


getFix.default <- function(fit, ...) stop(paste("Write a 'getFix' method for class",class(fit)))

#
#  getData is to lme what model.frame is to lm
#  getData is implemented as a method in nlme
#  so we just add a method for 'lm' objects
#

getData <- function(x,...) UseMethod("getData")
getData.lmer <- function(x,...) x@frame
getData.lme <- function(x,...) nlme:::getData.lme(x,...)
getData.lm <- function(x,...) model.frame(x,...)


# get the names of variables that are factors

getFactorNames <- function(object, ...) UseMethod("getFactorNames")

getFactorNames.data.frame <- function(object,...) {
     names(object)[ sapply(object, is.factor)  ]
}

getFactorNames.default <- function(object,...) getFactorNames( getData(object))


## print method for objects of class 'cat'

print.cat <- function(object,...) {
    cat(object)
    invisible(object)
}

# 
Lform <- function( fit, form, data = getData(fit)) {
 # 2011-12-01: replaced with version below
      if (missing(form)) return ( Lcall(fit))
      gg <- getFix(fit)
      Lsub <- do.call(cbind,eval( substitute( form ), data))
      if( (nrow(Lsub) != nrow( data))) {
          if ((nrow(Lsub)==1)) Lsub <- Lsub[rep(1,nrow(data)),]
          else stop('nrow(Lsub) != nrow(data)')
      }
      if( is.null( colnames(Lsub))) colnames(Lsub) <- rep('',ncol(Lsub))
      L <- matrix( 0, nrow = nrow(Lsub), ncol = length( gg$fixed))
      rownames(L) <- rownames(data)
      colnames(L) <- names( gg$fixed)
      Lpos <- Lsub[, colnames(Lsub) == '', drop = FALSE]
      # disp(Lpos)
      Lnamed <- Lsub[ , colnames(Lsub) !='', drop  = FALSE]
      # disp(Lnamed)
      for ( ip in seq_len( ncol( Lpos ))) L[,ip] <- Lpos[,ip]
      if ( ncol( Lnamed ) > 0 ) {
         if ( length( unknown <- setdiff( colnames(Lnamed) , colnames(L)))) {
              stop( paste("Unknown effect(s):" , unknown, collapse = " "))
         }
         for ( nn in colnames(Lnamed)) L[,nn] <- Lnamed[,nn]
      }
      attr(L,"data") <- data
      L
}

# 2012 12 04
# Plan for Lform
# 

# 2012 12 05: Lform becomes Lex to acknowledge the fact that it uses
# expressions instead of formulas



Lmat <- function(fit, pattern, fixed = FALSE, invert = FALSE, debug = FALSE) {
     # pattern can be a character used as a regular expression in grep
     # or a list with each component generating  a row of the matrix
     umatch <- function( pat, x ) {
            ret <- rep(0,length(pat))
            for ( ii in 1:length(pat)) {
                imatch <- grep(pat[ii], x, fixed= fixed, invert = invert)
                if ( length(imatch) != 1) {
                   cat("\nBad match of:\n")
                   print(pat)
                   cat("in:\n")
                   print(x)
                   stop("Bad match")
                }
                ret[ii] <- imatch
            }
            ret
     }
     if ( is.character(fit)) {
        x <- pattern
        pattern <- fit
        fit <- x
     }
     fe <- getFix(fit)$fixed
     ne <- names(fe)
     if (is.character(pattern)) {
        L.indices <- grep(pattern,names(fe))
        ret <- diag( length(fe)) [L.indices,,drop = FALSE]
        if (debug) disp(ret)
        rownames(ret) <- names(fe) [L.indices]
        labs(ret) <- "Coefficients"
     } else if (is.list(pattern)){
        ret <- matrix(0, nrow = length(pattern), ncol = length(fe))
        colnames(ret) <- ne
        for ( ii in 1:length(pattern)) {
            Lcoefs <- pattern[[ii]]
            pos <- umatch(names(Lcoefs), ne)
            if ( any( is.na(pos))) stop("Names of L coefs not matched in fixed effects")
            ret[ii, pos] <- Lcoefs
        }
        rownames(ret) <- names(pattern)
      }
      labs(ret) <- "Coefficients"
      ret
}

Ldiff.old <- function(fit, pat, levnames = c(reflevel,substring(rownames(L),cut+1)),
                       reflevel = "<ref>", cut = nchar(pat)) {
         L <- Lmat(fit, pat)
         nam <- rownames(L)
         n <- nrow(L)
         if(n < 2) return(L)
         plus <- unlist( apply( rbind( 2:n), 2, seq, n))
         minus <- rep(1:(n-1), (n-1):1)
         Lp <- L[ plus, ]
         Lm <- L[ minus, ]
         Lret <- rbind( L, Lp - Lm)
         rn <- paste( levnames[ c(1:n,plus) + 1], levnames[ c(rep(0,n), minus)+1], sep = " - ")
         rownames(Lret) <- rn
         Lret
}

Ldiff.rdc <- function( fit, nam , ref = "no longer used") {
      # based on Lmu
      # Tests all pairwise difference in factor with model with Intercept term
      Lm <- Lmu(fit, nam)
      levs <- rownames(Lm)
      n <- nrow(Lm)
      if (n < 2) return (Lm)
      plus <- unlist( apply ( rbind(2:n), 2, seq, n))
      minus <- rep(1:(n-1), (n-1):1)
      Lret <- Lm[plus,] - Lm[minus,]
      rn <- paste( levs [plus], levs[minus] , sep = " - ")
      rownames(Lret) <- rn
      Lret
}



Ldiff <- function( fit, pat, levnames = c(reflevel,substring(rownames(L),cut+1)),
         reflevel ="<ref>", cut=nchar(pat),verbose=F) {
      L <- Lmat(fit, paste("^",pat,sep=""))
      nam <- rownames(L)
      n <- nrow(L)
      zm <- matrix(1:n,nrow=n,ncol=n)
      plus <- zm[ col(zm) < row(zm)]
      minus <- rep(1:(n-1), (n-1):1)
      Lp <- L[plus,]
      Lm <- L[minus,]
      Lret <- rbind( L, Lp - Lm)
         pnames <- levnames [ c(1:n, plus) +1]
      mnames <- levnames [ c(rep(0,n), minus) + 1]
      if (verbose) {
        print(levnames)
        print(plus)
        print(minus)
        print(Lp)
        print(Lm)
        print(L)
        print(Lret)
        print(pnames)
        print(mnames)
      }
      rn <- paste( levnames[ c(1:n,plus)+1], levnames[ c(rep(0,n),minus) + 1], sep = " - ")
      rownames(Lret) <- rn
      Lret
}
 


Lmu <- function(fit, nam, verbose = 0) {
       ## "Works only if 'nam' is a factor and a main effect and model has Intercept")
       if ( class(fit) != 'lmer' ) stop( "only implemented for lmer")
       v <- fit@frame[[nam]]
       if( !is.factor(v)) stop ("nam needs to specify the name of a factor")
       levs <- levels(v)
       if( verbose > 0) print(levs)
       cmat <- contrasts(v)
       if( verbose > 0) print(cmat)
       #  print(cmat)
       fe <- getFix(fit)$fixed
       if( verbose > 0) print(fe)
       if ( substring(nam,1,1) != '^') nam <- paste("^",nam,sep="")
       L.indices <- grep(nam,names(fe))
       if( verbose > 0) print(L.indices)
       L <- matrix(0,nrow=length(levs), ncol = length(fe))

       colnames(L) <- names(fe)
       if( verbose > 0) print(L)
          rownames(L) <- levs
       L[,L.indices] <- cmat
       if('(Intercept)' %in% colnames(L)) L[,'(Intercept)'] <- 1
       L
}



Lc <- function(fit, nam, ref = 1, verbose = 0) {
       ## Comparisons with one level
       ## Use Lmu
       ## "Works only if 'nam' is a factor and a main effect and model has Intercept?")
       if ( class(fit) != 'lmer' ) stop( "only implemented for lmer")
       L <- Lmu( fit, nam)
       Lref <- L[ ref,,drop = FALSE]
       index <- 1:nrow(L)
       names(index) <- rownames(L)
       refind <- index[ref]
       if (length(refind) != 1) stop( paste( ref, "does not refer to a single level"))
       Lret <- L[-refind,]
       Lret <- Lret - cbind( rep(1,nrow(Lret))) %*% Lref
       attr(Lret,"heading") <- paste("Comparisons with reference level:", rownames(L)[refind])
       Lret
}

Lrm <- function(fit, nam, vals = 1:nrow(L.mu)) {
    ## Repeated measures polynomial contrasts
    ## Uses Lmu
    ## "Works only if 'nam' is a factor and a main effect and model has Intercept?")
    ##
    L.mu <- Lmu(fit, nam)
    # print(L.mu)
    pp <- cbind( 1, Poly(vals, nrow(L.mu) -1))
    # print(pp)
    ortho <- Q(pp)[,-1] # (linear, quad, etc.)
    # print(ortho)
    ortho <- ortho[,-1]
    maxp <- max( 5, nrow(L.mu))
    colnames(ortho) <- c('linear','quadratic','cubic', paste("^",4:maxp,sep=''))[1:ncol(ortho)]
    L <- t(ortho) %*% L.mu
    L
}
# Lrm(fit, "SRH94")

Lcall <- function( fit , factors = getFactorNames(fit), debug = F){

      nams <- names(getFix(fit)$fixed)

      nams <- gsub( "^", ":", nams)   # delineate terms
      nams <- gsub( "$", ":", nams)   # delineate terms
      for ( ff in factors)   {
          ff.string <- paste( ff, "([^:]*)" , sep = '')
          if(debug) disp( ff.string)
          ff.rep <- paste(ff, " == \\'\\1\\'", sep = '')
          if(debug) disp(ff.rep)
          nams <- gsub( ff.string, ff.rep, nams)
      }
     # for ( ii in seq_along(matrix)) {
     #     mm.all   <- paste( "(:",names(matrix)[ii], "[^\\)]*\\))",sep='')
     #     mm.match <- paste( "(",names(matrix)[ii], "[^\\)]*\\))",matrix[ii], sep ='')
     #     mm.rep   <- paste( "\\1")
     #     which.null <- grepl( mm.all, nams) mm.null  <-
     #
     # }
      nams <- sub("(Intercept)", 1, nams)
      nams <- gsub( "^:","(",nams)
      nams <- gsub( ":$",")",nams)
      nams <- gsub( ":", ") * (", nams)
      #if(comment) nams <- paste( nams, "  #",nams)
      nams <- paste( "with (data, \n cbind(", paste( nams, collapse = ",\n"), ")\n)\n", collapse = "")
      class(nams) <- 'cat'
      nams
}



Lequal <- function(fit, pat) {
       # Test for equality within levels of pat using all differences
         L <- Lmat(fit, pat)
         nam <- rownames(L)
         n <- nrow(L)
         if(n < 2) return(L)
         plus <- unlist( apply( rbind( 2:n), 2, seq, n))
         minus <- rep(1:(n-1), (n-1):1)
         Lp <- L[ plus, ]
         Lm <- L[ minus, ]
         Lret <- rbind( Lp - Lm)
         rn <- paste( nam[plus], nam[minus], sep = " - ")
         rownames(Lret) <- rn
         Lret
}




#Lc <- function(fit, vec ){
#   fe <- getFix(fit)$fixed
#   ret <- 0 * fe
#   if ( is.null(names(vec))) ret[] <-
#}
# Lmu(fit,"SRH")

Lall <- function( fit , nam ) {
        if ( class(fit) != 'lmer' ) stop( "only implemented for lmer")
        v <- fit@frame[[nam]]
        if( !is.factor(v)) stop ("nam needs to specify the name of a factor")
        lev0 <- levels(v)[1]
        ret <-list()
        namf <- nam
        if ( substring(namf,1,1) != "^") namf <- paste("^", namf, sep ="")
        ret[[ nam ]] <- Lmat( fit, namf)
        ret[[ paste(nam,"mu",sep = '.') ]] <- Lmu(fit,nam)
        ret[[ paste(nam,"diff",sep = '.') ]] <- Ldiff( fit , nam)
        ret

}


xlevels <- function(f,   type = c("raw","<mean>"),
        all = FALSE, sep = '') {
      # produces the equivalent of differentiation for a factor
      # optionally pairwise differences,
      # differences from the mean of other levels
      # or the difference from the weighted mean of other levels

      # BUG?: assumes no hyphens in original factor -> change sep

      if( !is.factor (f))   f <- factor( f, levels = f)
      cmat <- contrasts(f)
      nams <- levels(f)
      n <- length(nams)
      Pmats <- lapply(  type, function(x) Pmat( f, x, all = all, sep = sep))
      Pmat <- do.call( rbind, Pmats)
      #disp(Pmat)
      #disp(cmat)
      Cmat <- Pmat %*% cmat
      if( length(unique(rownames(Cmat))) == 1)  single <- TRUE else single <- FALSE
      if (single) Cmat <- rbind( Cmat,'___' = 0 )     #   CHANGE THIS
      fr <- factor(rownames(Cmat), levels=unique(rownames(Cmat)))
      contrasts(fr,n-1) <- Cmat
      names(fr) <- fr
      if (single) fr <- fr[-length(fr)]
      fr
}



dlevels <- function(f, type = "pairwise", all = FALSE, sep = '') {
      # produces the equivalent of differentiation for a factor
      # optionally pairwise differences,
      # differences from the mean of other levels
      # or the difference from the weighted mean of other levels

      # BUG?: assumes no hyphens in original factor -> change sep

      if( !is.factor (f))   f <- factor( f, levels = f)
      cmat <- contrasts(f)
      nams <- levels(f)
      n <- length(nams)
      Pmats <- lapply(  type, function(x) Pmat( f, x, all = all, sep = sep))
      Pmat <- do.call( rbind, Pmats)
      Cmat <- Pmat %*% cmat
      if( length(unique(rownames(Cmat))) == 1)  single <- TRUE else single <- FALSE
      if (single) Cmat <- rbind( Cmat,'<NULL>' = 0 )     #   CHANGE THIS
      fr <- factor(rownames(Cmat), levels=unique(rownames(Cmat)))
      contrasts(fr,n-1) <- Cmat
      names(fr) <- fr
      if (single) fr <- fr[-length(fr)]
      fr
}


# wtd mean of all but 1
# The common problem in all of the following is the generation of a Pmat matrix
# that combines rows of 'cmat' to achieve what it wants

Pmat <- function( f ,
      type = c("factor","raw","mean","(mean)","<mean>","II",
      "cen","cent","center","centre","<centre>" , "<center>" ,
      "(center)","(centre)","III","pairwise",  "<others.m>", "<others.c>",
        "diff", "diffmean","diffc","diffcen","diffcentre","diffcentre"),
        all = FALSE, sep = '')
        {
     if( length(type) > 1) return( do.call( rbind, lapply( type, function( t ) Pmat( f, t))))
     type <- match.arg(type)
     if( ! is.factor(f)) f <- factor( f, levels = unique(f))
     switch( type,
        raw =, factor = {
            ret <- diag(length(levels(f)))
            dimnames(ret) <- list(levels(f),levels(f))
            ret
        },
        cen =, cent =, center =, centre = , "III" = , "(centre)" = , "(center)", "<centre>" = , "<center>"  = {
            nlevs = length( levels(f) )
            matrix(rep(1/nlevs, nlevs), nrow = 1, dimnames = list(type, levels(f)))
        },
        mean =, "(mean)"=,  "<mean>"=, "II" = {
            matrix(table(f) / length(f), nrow = 1, dimnames = list(type, levels(f)))
        },
        pairwise = {
            nams <- levels(f)
            n <- length(nams)
            #if( length(nams) == 2) all = TRUE
            if ( all ) { # all comparisons
                plus <- rep( 1:n, n)
                minus <- rep (1:n, each = n)
                drop <- plus == minus
                plus <- plus[!drop]
                minus <- minus[!drop]
            } else {  # no duplicates
                zm <- matrix(1:n, nrow = n, ncol = n)
                plus <- zm[col(zm) < row(zm)]
                minus <- rep(1:(n - 1), (n - 1):1)
            }
            nrows <- length(plus)
            ret <- matrix( 0, nrows, n)
            ret[ cbind(1:nrows,minus) ] <- -1
            ret[ cbind(1:nrows,plus) ] <- 1
            rownames(ret) <- paste( nams[plus],'-', nams[minus], sep='')
            colnames(ret) <- nams
            ret
        } ,
        diff =, diffmean =, "<others.m>"= {
          ret <- Pmat_diffmat( table(f)/length(f) )
          rownames(ret) <- paste(rownames(ret),"-<others.m>", sep = sep)
          ret
        } ,
        diffcentre =, diffcenter =, diffc =, diffcen = , "<others.c>"= {
          levs <- levels(f)
          w <- rep(1/length(levs), length(levs))
          names(w) <- levs
          ret <- Pmat_diffmat( w )
          rownames(ret) <- paste(rownames(ret),"-<others.c>", sep = sep)
          ret
        })
}



Pmat_diffmat   <- function( w )  {
# generates matrix to apply to contrasts to yield differences
# from means of other levels (perhaps more interesting than pairwise?)
# A utitlity function for 'Pmat'
     ret <- - outer( 1/(1-w) , w, "*")
     diag(ret) <- 1
     dimnames( ret ) <- list( names(w), names(w))
     ret
}



wald.transform <- function (x, fun, label = "Transformed coefficients") {
  # transforms estimates and confidence intervals for wald, uses delta method
  # for se
  ant <- x[[1]]$estimate
  coefs <- as.double(ant$Estimate)
#   disp(coefs)
  trs <- numericDeriv( quote(fun(coefs)), 'coefs')
#   disp(trs)
  ant$Estimate <- c(trs)
  derivs <- abs( diag(as.matrix(attr(trs,'gradient'))))
  ant[["Std.Error"]] <- derivs * ant[["Std.Error"]]
  low.ind <- grep("^Lower ", colnames(ant))
  up.ind <- grep("^Upper ", colnames(ant))
  ant[[low.ind]] <- fun(ant[[low.ind]])
  ant[[up.ind]] <- fun(ant[[up.ind]])
  attr(ant,'labs') <- label
  x <- x[1]
  x[[1]]$estimate <- ant
  class(x) <- 'wald'
  x
}

