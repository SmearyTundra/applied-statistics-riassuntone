candle <- read.table("candle.txt", header = T)
sunshine <- read.table("sunshine.txt", header = T)
n <- 50
p <- 2
data <- rbind(candle, sunshine)
fact <- as.factor(c(rep("candle", 50), rep("sunshine", 50)))
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
library(mvnormtest)
mcshapiro.test(candle)
mcshapiro.test(sunshine)
S1 <- cov(candle)
S2 <- cov(sunshine)
S_pooled <- (49*S1 + 49*S2)/98
S1
S2
S_pooled
#we can assume gaussianity and homoschedasticity

smean_candle <- sapply(candle, mean)
smean_sunshine <- sapply(sunshine, mean)

#we perform a manova test
man <- manova(as.matrix(data) ~ fact)
summary.aov(man)
summary.manova(man)
#we have statistical evidence to reject the null hypothesis H0: mu_candle = mu_sunshine
mcshapiro.test(man$residuals)

#alternatively: T_stat
alpha <- 0.05
T_stat <- ((1/n + 1/n)^(-1))*(smean_candle - smean_sunshine)%*%solve(S_pooled)%*%(smean_candle - smean_sunshine)
cfr.fisher <- ((n + n - 2)*p/(n + n - 1 - p))*qf(1-alpha, p, n + n - 1 - p)
Rej <- T_stat > cfr.fisher
p_value <- 1 - df(T_stat/((n + n - 2)*p/(n + n - 1 - p)), p, n + n - 1 - p)

alphaB <- alpha/2
cfr.student <- qt(1- alphaB/2, n+n-2)

sdelta = smean_candle - smean_sunshine
IC <- cbind(inf = sdelta - cfr.student*sqrt(diag(S_pooled)*(1/n + 1/n)),
            center = sdelta,
            sup = sdelta + cfr.student*sqrt(diag(S_pooled)*(1/n + 1/n)))
rownames(IC) <- c("LM1 candle - sunshine", "LM2 candle - sunshine")
IC
#So with overall confidence 95% we can say that the candle lightbulbs are brighter on the short range, while the sunshine are brighter on the long range


#H0 : (LM1c - LM2c) - (LM1s - LM2s) > 0 vs H1
#we create a new dataset storing only the differences amongst the two measurements -> univariate
#

deltac <- candle[,1] - candle[,2] # ~ N(muc, a'SIGMAa)
deltas <- sunshine[,1] - sunshine[,2] # ~ N(mus, a'SIGMAa)
shapiro.test(deltac)
shapiro.test(deltas)
var.test(deltac, deltas)

#deltac ~ N(muc, sigma)
#deltas ~ N(mus, sigma)
#sdeltac - sdeltas ~ N(muc - mus, 2*sigma/n)
#sp = 49*s1 + 49*s2 / 98 ~ sigmaChi(98)/98
#(sdeltac - sdeltas)/sqrt(2*sp/n) ~ t(98)
x11()
boxplot(I(data[,1] - data[,2]) ~ fact)

sdeltac <- mean(deltac)
sdeltas <- mean(deltas)
sc <- sd(deltac)^2
ss <- sd(deltas)^2
sp <- ((n-1)*sc + (n-1)*ss)/(n+n-2)

T_stat <- (sdeltac - sdeltas)/sqrt(2*sp/n)
cfr.student.test <- qt(1-alpha, n+n-2)

Rej.test <- T_stat > cfr.student.test
p_val <- 1- pt(T_stat, n+n-2)


#alternatively
t.test(deltac, deltas, alternative= "greater", var.equal = T,
       conf.level = 1-alpha, paired = F, mu = 0 )

