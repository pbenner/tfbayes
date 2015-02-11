system("./distribution-test > distribution-test.csv")

t <- read.table("distribution-test.csv")$V1

hist(t,100,freq=F)
curve(dgamma(x, 0.5),0,8,add=T)
