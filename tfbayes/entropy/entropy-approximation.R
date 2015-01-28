
t2 <- read.table("entropy-approximation-2.csv")
t3 <- read.table("entropy-approximation-3.csv")
t4 <- read.table("entropy-approximation-4.csv")
t5 <- read.table("entropy-approximation-5.csv")

par(mfrow=c(2,2))
plot(t2, type="l", log="y", xlab="h", ylab="pdf")
plot(t3, type="l", log="y", xlab="h", ylab="pdf")
plot(t4, type="l", log="y", xlab="h", ylab="pdf")
plot(t5, type="l", log="y", xlab="h", ylab="pdf")
