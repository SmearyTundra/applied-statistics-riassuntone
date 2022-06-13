terrassa <- read.table("terrassa.txt", header = T)
girona <- read.table("girona.txt", header = T)
n <- 35
p<-2
#We are assuming that the datasets come from two bivariate gaussian distributions, we first
#need to verify wether we can assume the same covariance structure or not
St <- cov(terrassa)
Sg <- cov(girona)
St
Sg
var.test(terrassa$T2,girona$T2)
var.test(terrassa$T1,girona$T1)
#we can
#So we are performing the test on the difference of the mean of two bivariate gaussians
#it's a manova test
#H0: mu_terrassa = mu_girona vs H1
#We can hence perform a Fisher-distribution test
alpha <- 0.05
Sp <- ((n-1)*St + (n-1)*Sg)/(2*n-2)
cfr.fisher <- (p*(2*n - 2)/(2*n-1-p))*qf(1-alpha, p, 2*n-1-p)
smean.t <- sapply(terrassa, mean)
smean.g <- sapply(girona, mean)
T_stat <- (n/2)*(smean.t-smean.g)%*%solve(Sp)%*%(smean.t-smean.g)
Rej <- T_stat>cfr.fisher
Rej
#so we reject the null hypothesis, there is a significant difference in the mean, at level 95%

#Using Bonferroni correction we provide confidence interval for the difference in mean evaluation for each tester
#T.test <- a'(delta.smean) ~ N(a'delta.mean, a'(2/n SIGMA)a)
alphaB <- alpha/2
q.t <- qt(1-alphaB/2, 2*n - 2)
CI <- cbind(inf = smean.g - smean.t - q.t*sqrt(diag(Sp)*2/n),
            center = smean.g - smean.t,
            sup = smean.g - smean.t + q.t*sqrt(diag(Sp)*2/n))
rownames(CI) <- c("gvt_T1", "gvt_T2")
CI
#So girona has, on average, been better evaluated by both testers, in particular from the first one


#we are testing for a linear combination of the mean vectors:
#H0: mu.girona.T1 + mu.girona.T2 <= mu.terrassa.T1 + mu.terrassa.T2
#we can test it by using an unilateral t test
#T_test <- smean.gir.1 + smean.gir.2 - smean.gir.1 - smean.gir.2, reject when high

q.t.test <- qt(1-alpha, 2*n-2)
vec <- c(1,1)
T_test <- (vec%*%(smean.g - smean.t))/sqrt((2/n)*vec%*%Sp%*%vec)
Rej.test <- T_test > q.t.test
Rej.test
#yes, even if we discard who evaluated the tapas girona is still better evaluated than terrassa