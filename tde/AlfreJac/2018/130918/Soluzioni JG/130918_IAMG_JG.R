IAMG <- read.table("IAMG.txt", header = T)
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
mcshapiro.test(IAMG)
#we can assume a gaussian distribution
#The test statistic is:
#n*(smean - mean)'(S)^-1(smean - mean) ~ (n-1)*p/(n-p)*F(p, n-p)
#So the Ellipse is:
#{x e R3 s.t. (smean - x)'(S)^-1(smean - x) <= radius^2}
#radius = sqrt(((n-1)*p/(n*(n-p)))*quantile_fisher(1-alpha, p, n-p))
n <- dim(IAMG)[1]
p <- dim(IAMG)[2]
alpha <- 0.05


radius <- sqrt(((n-1)*p/(n*(n-p)))*qf(1-alpha, p, n-p))
smean <- sapply(IAMG, mean)
S <- cov(IAMG)

center <- smean
axes <- eigen(S)$vectors
lenghts <- sqrt(eigen(S)$values)*radius

#T2 confidence intervals are just the projections of this ellipse:
CI <- cbind(inf = smean - radius*sqrt(diag(S)),
            center = smean,
            sup = smean + radius*sqrt(diag(S)))
CI

#We want to test:
#H0: 0.1*mu_Registered - mu_No.show = 0 vs H1
a <- c(0.1, 0, -1)
#performing a t test

T_stat <- sqrt(n)*(abs(a%*%smean))/(sqrt(a%*%S%*%a))
cfr.student <- qt(1 - alpha/2, n-1)
Rej <- T_stat > cfr.student
p_value <- 1- pt(T_stat, n-1)
#and if the test were:
#H0: 0.1*mu_Registered - mu_No.show < 0 vs H1
T_stat_bis <- sqrt(n)*(a%*%smean)/(sqrt(a%*%S%*%a))
cfr.student_bis <- qt(1 - alpha, n-1)
Rej_bis <- T_stat_bis > cfr.student_bis
p_value_bis <- 1- pt(T_stat_bis, n-1)
