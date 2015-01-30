
source("entropy.R")
source("entropy-plot.R")

################################################################################

n  <- 100000
a1 <- 100
a2 <- 20

# execute sampler
system(sprintf("./entropy-distribution-test %d %d %d > entropy-distribution-test.csv", n, a1, a2))
# load result
samples <- read.table("entropy-distribution-test.csv")
samples.entropy <- apply(samples, 1, entropy)

simplex.contour(samples, filled=TRUE)

hist(samples.entropy, 100, prob=T, xlim=c(0, log(3)))
par(new=TRUE)
curve(dentropy(x, a1, a2), 0, log(3), xaxt="n", yaxt="n")
axis(4)
