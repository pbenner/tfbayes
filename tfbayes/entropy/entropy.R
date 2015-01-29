# basic definitions
################################################################################

entropy <- function(theta) -sum(sapply(theta, function(x) x*log(x)))
