
source("entropy.R")
source("entropy-plot.R")

################################################################################

dentropy <- function(h, k, a1 = 100, a2 = 20, ...) {
  dbeta(h/log(k), shape1 = a1, shape2 = a2)/log(3)
}

################################################################################

n  <- 100000
k  <- 3
a1 <- 100
a2 <- 20

# execute sampler
system(sprintf("./entropy-distribution-test %d %d %d %d > entropy-distribution-test.csv", n, k, a1, a2))
# load result
samples <- read.table("entropy-distribution-test.csv")
samples.entropy <- apply(samples, 1, entropy)

simplex.contour(samples, filled=TRUE)

hist(samples.entropy, 100, prob=T, xlim=c(0, log(k)))
par(new=TRUE)
curve(dentropy(x, k, a1, a2), 0, log(k), xaxt="n", yaxt="n")
axis(4)
