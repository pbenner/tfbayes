
library(KernSmooth)

# two-dimensional plotting functions
################################################################################

simplex.contour <- function(samples, bandwidth=c(0.02, 0.02), color=terrain.colors, filled=FALSE, ...) {
  if (dim(samples)[2] != 3) {
    return (NULL)
  }
  mycol <- rgb(0, 0, 0, alpha=0.5)
  est <- bkde2D(samples[,1:2], bandwidth=bandwidth, range.x=list(c(0,1),c(0,1)))
  if (filled) {
    with(est, filled.contour(x1, x2, fhat, color=color, ...))
  }
  else {
    with(est, contour(x1, x2, fhat, drawlabels=FALSE, col=mycol, ...))
    lines(c(0,1),c(1,0))
  }
}

simplex.plot <- function(f, n=100, col=rgb(0, 0, 0, alpha=0.5), color=terrain.colors, filled=FALSE, ...) {
  g <- function(x) {
    if (1-sum(x) <= 0) {
      return (0)
    }
    f(append(x, 1-sum(x)))
  }
  x  <- 1:n/n
  y  <- 1:n/n
  z  <- matrix(apply(expand.grid(x, y), 1, g), nrow=n, ncol=n)
  if (filled) {
    filled.contour(x, y, z, color=color, axes=FALSE, frame.plot=FALSE, ...)
  }
  else {
    contour(x, y, z, drawlabels=FALSE, col=col, axes=FALSE, frame.plot=FALSE, ...)
    lines(c(0,0,1,0),c(0,1,0,0))
  }
}
