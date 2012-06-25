
x1 <- c(1, 3, 5, 7, 9, 11)
y1 <- c(1, 5.243, 16.424, 40.182, 85.783, 163.919)

x2 <- c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39)
y2 <- c(1, 1.244, 1.697, 2.301, 3.16, 4.334, 5.816, 7.702, 10.286, 13.629, 17.311, 22.701, 28.121, 36.677, 46.119, 56.855, 70.592, 85.545, 104.637, 127.678)

xlim <- c(0, 40)
ylim <- c(0, max(y2)-10)

pdf("phylotree-expected-terms.pdf")

plot(xlim, ylim, type='n', xlab="Number of nodes", ylab="Expected number of terms")
lines(x1, y1, type='l')
lines(x2, y2, type='l')
text(12, 60, "Felsenstein's\nalgorithm")

dev.off()
