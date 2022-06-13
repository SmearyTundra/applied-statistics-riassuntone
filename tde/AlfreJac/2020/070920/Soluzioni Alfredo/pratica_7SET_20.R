## September 7th 2020

################
## PRoblem 2
#################
setwd("~/Documents/Polimi/APPSTATS")
candle <- read.table("candle.txt")
sunshine <- read.table("sunshine.txt")
head(candle)
head(sunshine)

# a)
# note samples are iid, and two distributions are independent
# we are working at n1 = 50, n2 = 50 with p = 2
"We could perform an asymptotic test by saying the numbers are large. 
We will be a little pignoli and say that we cannot apply LLN just yet"

"Hence, we will apply hotelling's theorem to see if these means differ.
We then assume independent multivar populations, with equal covariances"

# we first test normality
library(mvnormtest)
mshapiro.test(t(candle))
mshapiro.test(t(sunshine))
# with p values we do not reject null hypothesis

# now homoschedasticity: inspect visually covariances
S1 <- cov(candle)
S2 <- cov(sunshine)
x11(width=21)
par(mfrow=c(1,2))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
# to go on with analysis, we say they are approximately equal
# so no fisher-beherns problem

# Test
# H0: mu2 - mu1 = 0 H1 mu2 - mu1 != 0
#
delta0 <- rep(0, dim(candle)[2])
# Estimate for mus are the sample means, each distributed as Np(mu, 1/nSigma)
# with n the corresponding number of samples for each group
# so their difference will have mean mu2 - mu1, covariance (1/n1 + 1/n2) Sigma
xbar1 <- sapply(candle, mean)
xbar2 <- sapply(sunshine, mean)
sdiff <- xbar2 - xbar1
# To estimate the covariance, we get the pooled sample covariances estimate,
# which is distributed as Wishart( sigma / n1+n2-2, n1 +n2 - 2 df)

n1 <- dim(candle)[1]
n2 <- dim(sunshine)[1]
p <- dim(candle)[2]
Sp <- ((n1 -1 ) * cov(candle) + (n2-1) * cov(sunshine))/(n1+n2-2)
Spinv <- solve(Sp)
# Knowing both those distributions, we apply Hotelling's theorem
# We build the T2 statistic
T2 <- (1/n1 + 1/n2)**-1 * t(sdiff - delta0) %*% Spinv %*% ((sdiff - delta0))
# which is distributed, according to H's theorem,
# as (n1+n2-2)p/(n1+n2-1-p) F(p, n1+n2-1-p)

# therefore we reject at confidence level 95% if under the null hypothesis 
# that that the px1 0 vector is the true mean by seeing how likely is the T2 value
# assuming the null hypothesis

# confer the soglia from the F distribution, after which T² is too big to ignore
cfr.F <- p*(n1+n2-2)/(n1+n2-1-p)* qf(0.95, p, n1+n2-1-p)
T2 < cfr.F

# the T2 is very (less than the soglia) unlikely under the
#null hypothesis so we reject the null hypothesis

# b)
# Compute the p value:
pval <- 1 - pf(T2 / ( p*(n1+n2-2)/(n1+n2-1-p) ), p, n1+n2-1-p) 
pval
# we see indeed the probability under null hypothesis this value appears is lower
# than the soglia we fix at the requested level for the test

# c) 
# Produce Bonferroni confidence intervals for each marginal mean diff
C <- matrix(c(1, 0, 
              0, 1), 2, 2, byrow=T) 
k <- nrow(C)
# get 'radius' of t distribution
alpha <- .05
cfr.t <- qt(1-alpha/2/k ,n1 + n2 - 2)
# for first dimension differece
CI1 <- cbind(C[1,] %*% sdiff - cfr.t * sqrt(t(C[1,]) %*% Sp %*% C[1,] * (1/n1 + 1/n2)),
C[1,] %*% sdiff,
C[1,] %*% sdiff + cfr.t * sqrt( t(C[1,]) %*% Sp %*% C[1,]* (1/n1 + 1/n2)) )
CI2 <- cbind(C[2,] %*% sdiff - cfr.t * sqrt(t(C[2,]) %*% Sp %*% C[2,] * (1/n1 + 1/n2)),
             C[2,] %*% sdiff,
             C[2,] %*% sdiff + cfr.t * sqrt( t(C[2,]) %*% Sp %*% C[2,]* (1/n1 + 1/n2)) )
CI1
CI2
" Both confidence intervals contain the mu with 95% confidence, and since they don't cover
0, we reject the null hypothesis that the mu is 0, taking into
account we are performing simultaneous CIs (mantaining global alpha)"

# d) Compute the decreases for each brand, see if for Candle it is in mean higher
# than Sunshine
# H0: mu2 - mu1 = 0 H1 mu2 - mu1 < 0
# Compare wrt the t distribution but 1 sided: 
Dcandle <- candle[,1] - candle[,2]
Dssh <- sunshine[,1] - sunshine[,2]
bartlett.test(rbind(Dcandle, Dssh), rbind(rep("1", n1), rep("2", n2)))
# variances are the same too, as a corollary of the previous check 
t.test(Dcandle, Dssh, alternative= "greater", var.equal = T,
            conf.level = .95, paired = F, mu = 0 )
# param greater means candle has higher mean than sunshine
"We have statistical evidence to reject the null hypothesis in favor of the 
alternative which is that the decrease in candle is in mean higher than that
of sunshine, at confidence level 95%"
#############
# BONUS: approfondimento
#############
# do it manually to validate and learn
m1 <- mean(Dcandle)
m2 <- mean(Dssh)
sp2 <- ((n1 -1) * var(Dcandle) + (n2-1) * var(Dssh)) / (n1+n2-2)


"Parenthesis: a thought. Here we have a linear combination a = (1, -1)',
so we could get the same variance by applying the theorem of linar combinations
of gaussians to get each sd from its covariance. IT would be interesting to compute
the formula to obtain Sp in a similar way" # TODO

# confer the t distribution (one sided, higher tail)
cfr.t <- qt(0.95, df = n1 + n2 - 2)
tstatistic <- (m1 - m2)/sqrt(sp2 * (1/n1 + 1/n2))
tstatistic
pval <- 1 - pt(statistic, n1 + n2 -2)
pval
# CI
CIdiff <- cbind(inf = m1 - m2 - cfr.t * sqrt(sp2) * sqrt(1/n1 + 1/n2),
                center = m1 - m2,
                sup = m1 - m2 + cfr.t * sqrt(sp2) * sqrt(1/n1 + 1/n2 ))
CIdiff 
# Use hotelling's theorem (as T2 generalizes the univariate case)
# and the wishart generalizes chi quadro.
T2 <- (1/n1 + 1/n2)**-1 * (m1 -m2) * sp2**-1 * (m1-m2)
# First gauge check: t is sqrt(T²):
tstatistic
sqrt(T2)
T2
cfr.F <-  1*(n1+n2-2)/(n1+n2-1-1)* qf(0.95, 1, n1+n2-1-1)
cfr.F
1 - pf(T2, 1, n1 + n2 - 1-1) # we also reject: T2 makes the ellipse so it
# includes the worst direction onto which to project the data

"NB: it is possible to do simultaneous CI making it a one-sided test."