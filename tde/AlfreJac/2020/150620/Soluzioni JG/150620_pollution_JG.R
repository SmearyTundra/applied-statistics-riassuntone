pollution <- read.table("pollution.txt", header = T)
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
mcshapiro.test(pollution)
#we can assume gaussianity
#So ve are performing a test on the mean of a gaussian distribution:
#H0: mu = mu_0  vs H1
#n*(smean - mu_0)*S^-1*(smean - mu_0) ~ (n-1)*p/(n-p)*F(p, n - p + 1)
n <- dim(pollution)[1]
p <- dim(pollution)[2]
alpha <- 0.05

cfr.fisher <- ((n-1)*p/(n-p))*qf(1 - alpha, p, n - p)
smean <- sapply(pollution,mean)
mu_0 <- c(50,50)
S <- cov(pollution)

T_stat <-n*(smean - mu_0)%*%solve(S)%*%(smean - mu_0)

Rej <- T_stat > cfr.fisher
Rej

p_val <- 1 - pf(T_stat/((n-1)*p/(n-p)),p,n-p)
p_val

#ellipse
#{mu e R2 : (smean - mu)*S^-1*(smean - mu) <= radius^2}
radius <- sqrt(cfr.fisher/n)
directions <- eigen(S)$vectors
axes <- 2*sqrt(eigen(S)$values)*radius
library(car)
x11()
plot(pollution)
ellipse(smean, S, radius, col = 'blue')
x11()
plot(pollution, xlim = c(100, 200), ylim = c(100, 200))
ellipse(smean, S, radius, col = 'blue')

#mu_0 is heavily outside the confidence region
#indeed the previous test had strong evidence toward rejecting H0


CI <- cbind(inf = smean - radius*sqrt(diag(S)),
            center = smean,
            sup = smean + radius*sqrt(diag(S)))
CI