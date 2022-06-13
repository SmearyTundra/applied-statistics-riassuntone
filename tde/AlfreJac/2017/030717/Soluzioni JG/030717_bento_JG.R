bento <- bento[,-1]#remove the tags
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
library(mvnormtest)

#test for Gaussianity
mcshapiro.test(bento)
#not the best, but we can assume gaussianity of our data

#we want to perform the test: H0: mu(hanami) - mu(nohanami) = 0
#it's a test on a linear combination of the mean for a gaussian distribution
#we use matrix C =
#[1 0 0 0 -1 0 0 0]           [delta_rice]
#[0 1 0 0 0 -1 0 0]           [delta_sashimi]
#[0 0 1 0 0 0 -1 0] -> Cmu =  [delta_vegetables]
#[0 0 0 1 0 0 0 -1]           [delta_okashi]

#test stats: a'mu - 0/(sqrt(a'Sa/n)) ~ t(n-1)
#we use bonferroni, to test all 4 combinations simultaneously

smean <- colMeans(bento)
svar <- cov(bento)
n <- length(bento$rice_hanami)
a1 <- c(1, 0, 0, 0, -1, 0, 0, 0)
a2 <- c(0, 1, 0, 0, 0, -1, 0, 0)
a3 <- c(0, 0, 1, 0, 0, 0, -1, 0)
a4 <- c(0, 0, 0, 1, 0, 0, 0, -1)

T_stat_1 <- abs(a1%*%as.vector(smean)/(sqrt(t(a1)%*%svar%*%a1/n)))
T_stat_2 <- abs(a2%*%as.vector(smean)/(sqrt(t(a2)%*%svar%*%a2/n)))
T_stat_3 <- abs(a3%*%as.vector(smean)/(sqrt(t(a3)%*%svar%*%a3/n)))
T_stat_4 <- abs(a4%*%as.vector(smean)/(sqrt(t(a4)%*%svar%*%a4/n)))

alpha_Bonf <- 0.05/4
T_ref <-qt(1 - alpha_Bonf, n-1)
#test: we reject H0 if at least one of the T_stat > T_ref
(T_stat_1 > T_ref || T_stat_2 > T_ref || T_stat_3 > T_ref || T_stat_4 > T_ref)
#we reject H0 at level alpha


#we also verify that there is no strong evidence to assume correlation
x11()
image(svar)
#as expected, since the data collector provided independence


#SimCI(a'mu) = [a'X_bar +- sqrt(((n-1)*p/(n-p))*Quantile_Fisher_1-alpha(p,n-p))*sqrt(a'Sa/n)]
alpha <- 0.05
p <- dim(bento)[2]
S_inv <- solve(svar)


radius <- sqrt(((n-1)*p/(n-p))*qf(1-alpha, p, n-p))

C <- t(cbind(a1,a2,a3,a4))
T2 <- cbind(inf = C%*%as.vector(smean) - radius*diag(sqrt(C%*%svar%*%t(C)/n)),
            center = C%*%as.vector(smean), 
            sup = C%*%as.vector(smean) + radius*diag(sqrt(C%*%svar%*%t(C)/n)))
#so we notice that we do not evidence to support that there is a difference in the consumption of rice and vegetables
#however we have a strong evidence to support the fact that under Hanami families consume more sashimi and okashi
