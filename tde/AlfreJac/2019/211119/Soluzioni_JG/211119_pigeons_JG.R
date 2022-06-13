pigeons <- read.table("pigeons.txt", header = T)
data <- pigeons[, 1:2] - pigeons[, 3:4]
colnames(data) <- c("delta_weight", "delta_wing")
attach(data)
n <- dim(data)[1]
p <- dim(data)[2]

#we check if the data are gaussian
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
mcshapiro.test(data)
#we can assume it

#H0: mu_delta_weight = 0 & mu_delta_wing = 0 vs H1
smean <- colMeans(data)
S <- cov(data)

T_stat <- n*(smean)%*%solve(S)%*%smean/(n*(p-1)/(n-p))
p_val <- 1 - pf(T_stat, p, n-p)
p_val
#so we hace strong evidence to reject at any level
#indeed if we look at the Bonferroni corrected confidence intervals
alphaB <- 0.1/2
quant.t <- qt(1-alphaB/2,n-1)
CI <- cbind(inf = smean - quant.t*sqrt(diag(S)/n),
            center = smean,
            sup = smean + quant.t*sqrt(diag(S)/n))
CI

#if we were to interpret just the two populations:
males <- pigeons[,1:2]
females <- pigeons[,3:4]

mcshapiro.test(males)
mcshapiro.test(females)

S1 <- cov(males)
S2 <- cov(females)
S_p <- ((n-1)*S1 + (n-1)*S2)/(n*2 - 2)
smeanm <- colMeans(males)
smeanf <- colMeans(females)
T_f <- (n/2)*(smeanm - smeanf)%*%solve(S_p)%*%(smeanm - smeanf)/((2*n-2)*p/(2*n-p-1))
p_val_mf <- 1 - pf(T_f, p, 2*n - 1 - p)
p_val_mf

#wings ratio
wings_ratio <- pigeons[,2]/pigeons[,4]
head(wings_ratio)
shapiro.test(wings_ratio)
x11()
qqnorm(wings_ratio)
qqline(wings_ratio)
#we can work with a gaussian approximation
#moreover we could even work asymptotically
#H0: wings_M/wings_F <= 1.2 vs H1
t.test(wings_ratio, mu = 1.2, alternative = "g")
#we have strong evidence supporting the claim that, on average, the wing span of a male pigeon is 1.20 the wing span of his companion