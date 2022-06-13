castle <- read.table("castle.txt", header = T)
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
library(mvnormtest)
library(car)

mcshapiro.test(castle)
#we can assume gaussianity
n <- dim(castle)[1]
p <- dim(castle)[2]

#Stat test : n*(sample_mean - mean)'(S)^-1(sample_mean - mean) ~ (n-1)*p/(n-p)*F(p, n-p)
#Test H0 : mean = mean0 vs H1
smean <- colMeans(castle)
S <- cov(castle)
mean0 <- c(45.733,7.333)
T_stat <- n*(smean - mean0)%*%solve(S)%*%matrix(smean - mean0)/((n-1)*p/(n-p))
p_val<- 1 - pf(T_stat, p, n-p)
p_val
#we cannot refute H0

#WE ARE NOT BUILDING A CONFIDENCE REGION FOR THE MEAN, BUT A PREDICTION REGION FOR THE OBSERVATION!!!
#x ~ N(mean, SIGMA) -> (X - s_mean) ~ N(0, (1+1/n)SIGMA) (sum of two gaussians) -> (X - s_mean)'inv((1+1/n)SIGMA)(X - s_mean) ~ Chisq(2) (but we do not know SIGMA)
#(X - s_mean)/sqrt((1+1/n)) ~ N(0, SIGMA), S ~ W(SIGMA/(n-1), n-1)
#((n-1)-p+1)/((n-1)*p)((X - s_mean)/sqrt((1+1/n)))'inv(S)((X - s_mean)/sqrt((1+1/n))) ~ F(p, (n-1) - p +1 )
#region {((n-1)-p+1)/((n-1)*p)((X - s_mean)/sqrt((1+1/n)))'inv(S)((X - s_mean)/sqrt((1+1/n))) <= c^2}
# = {(((n-1)-p+1)/(((n-1)*p)*(1+1/n)))*(X - s_mean)'inv(S)(X - s_mean) <= c^2}
alpha <- 0.05
radius <- sqrt(qf(1-alpha,p,n-p)/((n-p)/((n*p-p)*(1+1/n))))
center <- smean

directions <- eigen(S)$vectors
length <- eigen(S)$values*radius
#if we work asymptotically
radius_chi <- sqrt(qchisq(1-alpha, 2)*(1+1/n))

x11()
plot(castle, asp=1, pch=1, main='Dataset of the Differences')
ellipse(center=smean, shape=S,radius, lwd=2, col = 'blue')
ellipse(center=smean, shape=S,radius_chi, lwd=2, col = 'red')  
#similar results with using exact of asymptotic ellipse