
source("entropy.R")
source("entropy-plot.R")

################################################################################

dentropy <- function(h, a1 = 100, a2 = 20, ...) {
  dbeta(h/log(3), shape1 = a1, shape2 = a2)/log(3)
}

################################################################################

n  <- 100000
a1 <- 100
a2 <- 20

# execute sampler
system(sprintf("./entropy-distribution-test %d 2 %d %d > entropy-distribution-test.csv", n, a1, a2))
# load result
samples <- read.table("entropy-distribution-test.csv")
samples.entropy <- apply(samples, 1, entropy)

simplex.contour(samples, filled=TRUE)

hist(samples.entropy, 100, prob=T, xlim=c(0, log(3)))
par(new=TRUE)
curve(dentropy(x, a1, a2), 0, log(3), xaxt="n", yaxt="n")
axis(4)
