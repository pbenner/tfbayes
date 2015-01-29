
source("entropy-plot.R")

################################################################################

n  <- 100000
a1 <- 100
a2 <- 20

# execute sampler
system(sprintf("./entropy-distribution-test %d %d %d > entropy-distribution-test.csv", n, a1, a2))
# load result
theta <- read.table("entropy-distribution-test.csv")

simplex.contour(theta, filled=TRUE)
