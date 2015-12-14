# revised by MF to replace Splus calls 7/19/2010
# (removed brace from fun.R)





#' Create a brace for a graph
#' 
#' The function \code{brace()} calculates the (X,Y) coordinates to draw a brace
#' between two points in a graph.
#' 
#' 
#' @param x1,y1 Coordinate of the initial point for the brace
#' @param x2,y2 Coordinate of the final point for the brace
#' @param right Logical to indicate if a vertical brace should have its'
#' central point facing right. If the brace is horizontal (\code{y1 == y2}),
#' "right" means "down"
#' @param rad Radius of the quarter circles used to draw the brace
#' @return A 2-column matrix of (X,Y) points, suitable for use with
#' \code{lines}
#' @author Georges Monette, modified for R by Michael Friendly
#' @seealso \code{\link[graphics]{lines}}, \code{\link[graphics]{arrows}}, ~~~
#' @keywords dplot aplot
#' @examples
#' 
#' plot(c(-1,1), c(-1,1), type="n", xlab="x", ylab="y")
#' abline(h=0, col="gray")
#' abline(v=0, col="gray")
#' b <- 0.6
#' abline(0, b, col="blue")
#' lines(brace(0, 0, 0, b, right=FALSE, rad=0.1))
#' lines(brace(0, 0, 1, 0, rad=0.1))
#' text(0.5, -.2, '1', cex=2)
#' text(-.2, b/2, 'b', cex=2)
#' 
#' @export
brace <- function (x1 = 0, y1 = 0, x2 = 0, y2 = 1, right = TRUE, rad = 0.2) 
{
#   uin only in Splus
#    uin <- par("uin")
#    ux <- uin[1]
#    uy <- uin[2]
	pin <- par("pin")
	usr <- matrix(par("usr"), ncol=2)
	ux <- pin[1] / diff(usr[,1])
	uy <- pin[2] / diff(usr[,2])
	
	dx <- x2 - x1
	dy <- y2 - y1
#     atan(x,y) from Splus -> atan(a/y) in R 
#   alpha <- atan(ux * dx, uy * dy)
	alpha <- atan((ux * dx) / (uy * dy))
	scale <- sqrt((ux * dx)^2 + (uy * dy)^2)
	if (scale > 5 * rad) 
		rad <- rad/scale
	qcirc <- cbind(cos((0:10) * pi/20), sin((0:10) * pi/20))
	qcircr <- cbind(cos((10:0) * pi/20), sin((10:0) * pi/20))
	rot <- function(theta) t(cbind(c(cos(theta), sin(theta)), 
						c(-sin(theta), cos(theta))))
	seg1 <- t(t(rad * qcirc %*% rot(-pi/2)) + c(0, rad))
	seg4 <- t(t(rad * qcirc) + c(0, 1 - rad))
	seg3 <- t(t((rad * qcircr) %*% rot(pi)) + c(2 * rad, 0.5 + 
							rad))
	seg2 <- t(t((rad * qcircr) %*% rot(pi/2)) + c(2 * rad, 0.5 - 
							rad))
	bra <- rbind(seg1, seg2, seg3, seg4)
	if (!right) 
		bra <- bra %*% diag(c(-1, 1))
	bra <- scale * bra %*% rot(-alpha)
	bra <- bra %*% diag(c(1/ux, 1/uy))
	bra <- t(t(bra) + c(x1, y1))
	bra
}
