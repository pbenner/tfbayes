
source("entropy-plot.R")
source("entropy-density.R")

################################################################################

# counts
# ------------------------------------------------------------------------------

a1 <- 100
a2 <- 20

a1 <- 10
a2 <- 10

# sampler
# ------------------------------------------------------------------------------

samples <- rentropy.simplex(50000, 3, a1, a2, sampler.options=list(method = "dirichlet"))
samples.entropy <- apply(samples, 1, entropy)

samples <- rentropy.simplex(50000, 3, a1, a2, sampler.options=list(method = "normal", sd = 0.1))
samples.entropy <- apply(samples, 1, entropy)

samples <- rentropy.simplex(10000, 3, a1, a2, sampler.options=list(method = "normal"))
samples.entropy <- apply(samples, 1, entropy)

# plotting
# ------------------------------------------------------------------------------

hist(samples.entropy, 100, prob=T, xlim=c(0, log(3)))
par(new=TRUE)
curve(dentropy(x, a1, a2), 0, log(3), xaxt="n", yaxt="n")
axis(4)

simplex.contour(samples, filled=TRUE)

simplex.plot(function(x) dentropy.simplex(x, a1, a2), filled=TRUE)
