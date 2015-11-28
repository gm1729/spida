#' Lattice panel function for data ellipse
#' 
#' Draws a data ellipse in panels  
#' 
#' @param x 
#' @param y
#' @param ... 
#' 
#' @return called for its side effects of plotting
#' @examples
library(spida)
library(lattice)
library(latticeExtra)
# No groups or panels:
xyplot(mathach ~ ses, hs) + layer(gpanel.ellipse(...))
# Panels:
xyplot(mathach ~ ses | school, hs) + layer(gpanel.ellipse(...))
# Groups:
xyplot(mathach ~ ses, hs, groups = Sex) + glayer(gpanel.ellipse(...))
# Panels with groups:
xyplot(mathach ~ ses | school, hs, groups = Sex) + glayer(gpanel.ellipse(...))
#' @export
gpanel.ellipse <-
function(x, y, ..., type) {
  mat <- na.omit(cbind(x,y))
  if( nrow(mat) < 3) return(NULL) 
  ell <- dell(mat)
  panel.xyplot(ell[,1], ell[,2], ..., type = 'l')
}
