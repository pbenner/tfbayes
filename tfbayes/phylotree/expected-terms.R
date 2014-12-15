# required libraries
################################################################################
require(tikzDevice)

################################################################################

t <- read.table("expected-terms.dat", header=F)
colnames(t) <- c("x", "y1", "y2")
attach(t)

plot(x, y1, type='l',
     xlab="Number of vertices",
     ylab="Expected number of terms")
lines(x, y2, lty=2)
legend("bottomright", c("pruning algorithm", "likelihood decomposition"), lty=c(1,2), bty='n')

# export
################################################################################

tikz("expected-terms.tex", bareBones=TRUE, width=3.2, height=2.3)
par(mar=c(4.2,4,1,1))
plot(x, y1, type='l',
     xlab="Number of vertices",
     ylab="Expected number of terms")
lines(x, y2, lty=2)
legend("top", c("decomposition", "pruning"), lty=c(1,2), bty='n')
dev.off()
