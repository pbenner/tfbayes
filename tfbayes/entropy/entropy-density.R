
# approxfun
require(stats)
# kernel density estimates
require(ks)

source("entropy-sampler.R")

# compute volumes
################################################################################

n <- 1000000
k <- 3
volume.samples <- rdirichlet(n, rep(1, k))
volume.samples.entropy <- apply(volume.samples, 1, entropy)
rm(volume.samples)

volume.estimate <- kde(volume.samples.entropy, h=0.01, xmin=0, xmax=log(k))

volume <- function(theta) {
  predict(volume.estimate, x=entropy(theta))
}

## volume.histogram         <- hist(volume.samples.entropy, 100, plot=F)
## volume.histogram$mids    <- c(0, volume.histogram$mids,    log(k))
## volume.histogram$density <- c(0, volume.histogram$density, 0)

## volume <- function(theta) {
##   f <- approxfun(volume.histogram$mids, volume.histogram$density)
##   f(entropy(theta))
## }

# define densities
################################################################################

dentropy <- function(h, a1 = 100, a2 = 20, ...) {
  dbeta(h/log(3), shape1 = a1, shape2 = a2)/log(3)
}

dentropy.simplex <- function(theta, ...) {
  dentropy(entropy(theta), ...)/volume(theta)
}

rentropy.simplex <- function(n, k, ..., sampler.options=list()) {
  sampler.options$n <- n
  sampler.options$f <- function(x) dentropy.simplex(x, ...)
  sampler.options$k <- k
  do.call("simplex.sampler", sampler.options)
}
