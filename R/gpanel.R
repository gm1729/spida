## 
## panel functions:
##
## gpanel.fit incorporating gpanel.band
## gpanel.text
##
## designed to work easily with 'layer' and 'glayer' in latticeExtra
## Added: 2014 09 18
## Author: GM
##

# gpanel.fit ####
#   Adds error bands. Uses two additional arguments: 
#   'lower' and 'upper' for the lower and upper values of the bands
#   There are, provisionally, two versions:
#   panel.band which is minimal. To get group colors, specify col = col.symbol
#   and panel.band2 based on panel.smoother in latticeExtra that does a lot
#   more work with colors, etc. 
#   We'll see which approach works best.
#
#FUTURE: gpanel.model 
# gpanel.model( x, y, model, se = T)
# will combine gpanel.fit and gpanel.band
#


gpanel.fit <-
  function(x, y, fit, lower, upper,
           subscripts, ..., type, group.number, alpha, col, col.line, col.symbol, border = F, font, fontface)
  {
    if( !missing(fit)) {
      if( missing(col) ) col <- 'blue'
      if( !missing(group.number)) {
        col <- col.line
      }
      if( !missing(subscripts) ) {
        fit <- fit[subscripts]
      }
      dd <- data.frame( x=x, fit = fit)
      dd <- dd[order(dd$x),]
      panel.lines( dd$x, dd$fit, ..., col = col, type = 'l')      
    }
    if( !missing(lower)){
      if( missing(alpha) || alpha == 1) alpha <- .3
      if( missing(col) ) col <- 'blue'
      if( !missing(group.number)) {
        col <- col.symbol
      }
      if( !missing(subscripts) ) {
        upper <- upper[subscripts]
        lower <- lower[subscripts]
      }
      dd <- data.frame( x=x, lower = lower, upper = upper)
      dd <- dd[order(dd$x),]
      panel.polygon( c(dd$x, rev(dd$x)),c(dd$upper, rev(dd$lower)), 
                     border = border, col = col, alpha = alpha,...)
    }
      #  panel.polygon(c(dd$x, rev(dd$x)), c(dd$upper, rev(dd$lower)), col = col, alpha = alpha, ...)
  }


gpanel.band <- gpanel.fit
gpanel.labels <-
  function (x, y, labels , subscripts, ...) 
  {
    # NOTE: Also include names of anything you DON'T want passed.
    # this is an experiment in writing a function that can be
    # called via layer or glayer without further complications
    # e.g. xyplot(......,labels = rownames(data)) + layer( gpanel.labels(...))
    # or  xyplot(....., labels = rownames(data), subscripts = T) +
    # glayer(gpanel.labels(...))
    # for selected lables, use
    # trellis.focus()
    # panel.identify(labels = rownames(data),rot=-15,col = col.symbol, etc.)  
    #
    # Note that subscripting is taken care of
    # trellis.unfocus()
    
    if(!missing(subscripts)) {
      labels <- labels[subscripts]
    }
    panel.text(x, y, labels, ...)
  }


fillin <- function(data, form, n = 200, xpd = 1.0) {
  # data: data frame with values of x to fill in within ranges of x within
  #       levels of g
  # form: formula idenfying variable x to fill in and grouping variables, g1, g2, etc.
  #       of form ~ x + g1 + g2 (the variable to fill in comes first)
  # n:    number of values to fill in across the entire range
  # xpd:  proportional range expansion
  
  sel.mf <- model.frame(form, data, na.action = na.pass)
  ret <- sel.mf
  xmina <- min(sel.mf[[1]], na.rm = T)
  xmaxa <- max(sel.mf[[1]], na.rm = T)
  xmida <- (xmaxa+xmina)/2
  xrana <- (xmaxa-xmina)/2
  xmina <- xmida - xpd*xrana
  xmaxa <- xmida + xpd*xrana
  if (ncol(sel.mf) ==1 ) {
    ret$xmin <- min(sel.mf[[1]],na.rm=T)
    ret$xmax <- max(sel.mf[[1]],na.rm=T)
  } else {
    grps <- lapply( sel.mf[2:ncol(sel.mf)], function(x)
      as.numeric(as.factor(x)))
    grps <- factor( do.call(paste,grps))
    ret$xmin <- capply( sel.mf[[1]], grps, min, na.rm=T)
    ret$xmax <- capply( sel.mf[[1]], grps, max, na.rm=T)
  }
  summ <- up(ret, formula(ret)[-2])  # has every unique combination of g's
  added.vals <- summ[ rep(1:nrow(summ), n),]
  added.vals[[names(ret)[1]]] <- rep( seq(xmina, xmaxa, length.out = n), each = nrow(summ))
  x <- added.vals[[names(ret)[1]]]
  x.min <- added.vals[['xmin']]
  x.max <- added.vals[['xmax']]
  x.mid <- (x.min + x.max)/2
  x.ran <- (x.max - x.min)/2
  ok <- (x <= (x.mid + xpd*x.ran)) & (x >= (x.mid - xpd*x.ran))
  merge( data, added.vals[ok,], all = T)  
}

# zd <- data.frame( x = 1:10, 
#                   g = rep(c('a','b'), each = 5), 
#                   y = 1:10,
#                   g2 = LETTERS[c(1,1,1,2,2,2,3,3,3,3)])
# zd
# (ff <- fillin(zd, ~ x + g, n = 10, xpd = 1.2))
# tab(ff, ~ g)
# barchart
# ?barchart
# library(latticeExtra)
# zz <- fillin(zd, ~ x + g+g2, n =100)
# xyplot( x ~ g|g2 ,zd, cex = 2) +
#   xyplot( x ~ g|g2, fillin(zd, ~ x + g+g2, n =100),col='red')
# 
# as.data.frame(tab(zd, ~g))
