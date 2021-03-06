
source("entropy-plot.R")
source("entropy-density.R")

################################################################################

# counts
# ------------------------------------------------------------------------------

k <- 3

a1 <- 100
a2 <- 20

a1 <- 10
a2 <- 10

a1 <- 1
a2 <- 1

# visualize density
# ------------------------------------------------------------------------------

simplex.plot(function(x) dentropy.simplex(x, a1, a2), filled=TRUE)

# sampler
# ------------------------------------------------------------------------------

samples <- rentropy.simplex(50000, k, a1, a2, sampler.options=list(method = "dirichlet"))
samples.entropy <- apply(samples, 1, entropy)

samples <- rentropy.simplex(50000, k, a1, a2, sampler.options=list(method = "normal", sd = 0.1))
samples.entropy <- apply(samples, 1, entropy)

samples <- rentropy.simplex(10000, k, a1, a2, sampler.options=list(method = "normal"))
samples.entropy <- apply(samples, 1, entropy)

samples <- rentropy.simplex(10000, k, a1, a2, sampler.options=list(method = "hamilton", epsilon=0.0001, L=200, burnin=0))
samples.entropy <- apply(samples, 1, entropy)

samples <- rentropy.simplex(10000, k, a1, a2, sampler.options=list(method = "hamilton", epsilon=0.0005, L=50, burnin=0, permute=TRUE))
samples.entropy <- apply(samples, 1, entropy)

# plotting
# ------------------------------------------------------------------------------

hist(samples.entropy, 100, prob=T, xlim=c(0, log(k)))
par(new=TRUE)
curve(dentropy(x, k, a1, a2), 0, log(k), xaxt="n", yaxt="n")
axis(4)

simplex.contour(samples, filled=TRUE)

plot(samples[,1:2], type="l", xlim=c(0,1), ylim=c(0,1), ann=FALSE)
