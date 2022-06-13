olives <- read.table("olives.txt", header = T)
attach(olives)
n <- dim(olives)[1]
p <- 3
olives_L <- olives[which(Restaurant=="Dalla Luigina"),]
olives_M <- olives[which(Restaurant=="Caffè Muletti"),]
n_L <- dim(olives_L)[1]
n_M <- dim(olives_M)[1]

#we first of all need to check gaussianity
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
library(mvnormtest)
mcshapiro.test(olives_L[,1:3])
mcshapiro.test(olives_M[,1:3])
#we have gaussianity
S_L <- cov(olives_L[,1:3])
S_M <- cov(olives_M[,1:3])
S_L
S_M
x11()
image(1:p,1:p, S_L)
x11()
image(1:p,1:p, S_M)
bartlett.test(olives[,1:3], Restaurant)
#we can work with homoschedasticity

#the test we need to perform is a manova:
alpha <- 0.05
man <- manova(as.matrix(olives[,1:3])~olives[,4])
summary.manova(man)
summary.aov(man)
#we have strong evidence toward saying that there is a difference in the means of the two olives served
summary.manova(man)$stats[1,6] < alpha

#test statistic: (1/n1 + 1/n2)^(-1)[(smean1 - smean2) - delta0]'(S_pooled)^-1[(smean1 - smean2) - delta0] ~ ((n1 + n2 -2)*p/(n1 + n2 - 1 - p))*F(p, n1 + n2 - p)
#T2 intervals are just the projection of the overall confidence region
#SimCI(a'X) = [a'sample_delta_X +- sqrt(((n-1)*p/(n-p))*Quantile_Fisher_1-alpha(p,n-p))*sqrt(a'S_p a*(1/n1 + 1/n2))]
smean_L <- sapply(olives_L[,1:3], mean)
smean_M <- sapply(olives_M[,1:3], mean)
sdelta_L_M <- smean_L - smean_M 
S_pooled <- ((n_L - 1)*S_L + (n_M - 1)*S_M)/(n_L + n_M - 2)
cfr.fisher <- ((n_L + n_M -2)*p/(n_L + n_M - 1 - p))*qf(1-alpha,p, n_L + n_M - 1 - p)

increment <- sqrt(diag(S_pooled)*(1/n_M + 1/n_L))
T2CI <- cbind(inf= sdelta_L_M - cfr.fisher*increment, point = sdelta_L_M, sup = sdelta_L_M + cfr.fisher*increment)
rownames(T2CI) <- c("Tot_L - Tot_M", "Fill_L - Fill_M", "Bread_L - Bread_M")
T2CI
#Muletti makes on average bigger olives, with more filling
#there is no strong evidence to support the claim that it uses more bred, though