
# approxfun
require(stats)
# kernel density estimates
require(ks)

source("entropy-sampler.R")

# compute volumes
################################################################################

volume.data <- list(
            NULL,
            read.table("entropy-approximation-2.csv"),
            read.table("entropy-approximation-3.csv"),
            read.table("entropy-approximation-4.csv"),
            read.table("entropy-approximation-5.csv"))

volume.f <- function(h, k) {
  t <- volume.data[[k]]
  n <- nrow(t)
  v.min <- 0
  v.max <- log(k)
  v.width <- (v.max-v.min)/n
  i <- floor((h-v.min)/v.width)+1
  i <- if(i == n) i-1 else i
  # interpolate the result
  x1 <- t[i,  1];
  y1 <- t[i,  2];
  y2 <- t[i+1,2];
  y1 + (y2 - y1)*(h - x1)/v.width;
}

# define densities
################################################################################

dentropy <- function(h, k, a1 = 1, a2 = 1, ...) {
  dbeta(h/log(k), shape1 = a1, shape2 = a2)
}

dentropy.simplex <- function(theta, ...) {
  h <- entropy(theta)
  k <- length(theta)
  dentropy(h, k, ...)/volume.f(h, k)
}

rentropy.simplex <- function(n, k, ..., sampler.options=list()) {
  sampler.options$n <- n
  sampler.options$f <- function(x) dentropy.simplex(x, ...)
  sampler.options$k <- k
  do.call("simplex.sampler", sampler.options)
}
