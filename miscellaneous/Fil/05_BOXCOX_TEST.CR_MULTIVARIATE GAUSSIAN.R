### TOPIC:
### Box-Cox transformations
### Tests and confidence regions for the mean of a multivariate Gaussian

library(car)
library(mvtnorm)


################# STAT 101 TEST ##################
# Univariate
t.test(d, mu = 0.2, alternative = "greater")








################## BOX-COX #######################
## Univariate Box-Cox transformation
lambda.x <- powerTransform(x) 
lambda.x
# lambda<1: observations <1 are "spread", observations >1 are "shrinked"

# Transformed sample with the optimal lambda (command bcPower of library car)
bc.x <- bcPower(x, lambda.x$lambda)      # it transforms the data of the first argument through the Box-Cox 
                                         # transformation with lambda given as second argument
hist(bc.x, col='grey', prob=T, main='Histogram of BC(X)', xlab='BC(x)')
shapiro.test(x)
shapiro.test(bc.x)

## Multivariate Box-Cox (p=2)
lambda <- powerTransform(cbind(x,y))
lambda
# Compute the transformed data with optimal lambda (of the bivariate transf.)
# (command bcPower)
BC.x <- bcPower(x, lambda$lambda[1])
BC.y <- bcPower(y, lambda$lambda[2])

# MC-Shapiro test (H0: multivariate gaussianity)
mcshapiro.test(cbind(x,y))$p
mcshapiro.test(cbind(BC.x,BC.y))$p











################## TEST AND CONFIDENCE REGION FOR THE MEAN OF A MULTIVARIATE GAUSSIAN #######################
# 1)  Formulate the test (and test the Gaussian assumption, if needed)
# 2)  Compute the test statistics 
# 3a) Having set the level of the test, verify whether the test statistics 
#     belongs to the region of rejection (i.e., if there is statistical  
#     evidence to reject H0)
# 3b) Compute the p-value of the test

####### Test (GAUSSIANITY ASSUMED)
### Test on the mean of level alpha=1%
### H0: mu == mu0 vs H1: mu != mu0
### with mu0=c(1,0)
###-----------------------------------
mcshapiro.test(x)
x.mean   <- sapply(x,mean)
x.cov    <- cov(x)
x.invcov <- solve(x.cov)
n <- nrow(x) 
p <- ncol(x)  # CHANGE TO DIMENSION OF MULTIVARIETY 

alpha <- 0.01
mu0 <- c(1,0)

# T2 Statistics
x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
# Radius of the ellipsoid
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
# Test: 
x.T2 < cfr.fisher   # no statistical evidence to reject H0 at level alpha
# Rejection region: {x.T2>cfr.fisher}
# (we reject for large values of the T2 statistics)

# Compute the p-value 
P <- 1-pf(x.T2*(n-p)/((n-1)*p), p, n-p)
P


####### Region of confidence/rejection 
###----------------------------------
# Rejection region (centered in mu0)
# { m \in R^2 s.t. n * (mu0-m)' %*% (x.cov)^-1 %*% (mu0-m) < cfr.fisher }
x11()
plot(x, asp = 1)
ellipse(mu0, shape=x.cov/n, sqrt(cfr.fisher), col = 'blue', lty = 2, center.pch = 16)
points(x.mean[1], x.mean[2], pch = 16, col ='red', cex = 1.5) # We add a red point in correspondence of the sample mean

# Confidence region (centered in x.mean)
# { m \in R^2 s.t. n * (x.mean-m)' %*% (x.cov)^-1 %*% (x.mean-m) < cfr.fisher }
ellipse(x.mean, x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 2, lwd=2, center.cex=1)

# Note: by definition, the confidence region of level 1-alpha
# produces ellipsoidal regions that contain the true mean
# 100(1-alpha)% of the times if H0 is TRUE.

# Characterize region
# Center:
x.mean

# Directions of the principal axes:
eigen(x.cov/n)$vectors

# Length of the semi-axes of the ellipse:
r <- sqrt(cfr.fisher)
r*sqrt(eigen(x.cov/n)$values) 
# Warning: Conf Reg => x.cov/n   IF 99% CONF. REGION FOR DATA NOT MEAN x.cov WITHOUT /n





####### SIMULTANEOUS T2 CONFIDENCE INTERVALS
mu0   <- c(1,0)
alpha <- 0.01
x.mean   <- sapply(x,mean)
x.cov    <- cov(x)
x.invcov <- solve(x.cov)
n <- nrow(x) 
p <- ncol(x)  # CHANGE TO DIMENSION OF MULTIVARIETY 
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2

# Plot of confidence and rejection regions + T2 rectangle intervals
x11()
par(mfrow=c(1,1))
plot(x, asp = 1,main='Confidence and rejection regions')
ellipse(mu0, shape=x.cov/n, sqrt(cfr.fisher), col = 'blue', lty = 2, center.pch = 16)
points(x.mean[1], x.mean[2], pch = 16, col = 'red', cex=1.5)

ellipse(x.mean, shape=x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 2, center.pch = 16)
rect(T2[1,1],T2[2,1],T2[1,3],T2[2,3], border='red', lwd=2)


# Both the intervals contain the mean under H0
# (i.e., mu0 is contained in the rectangular region determined by
# the projection of the ellipsoid along the coordinate directions)
# Remark: this is not in contrast with the previous findings
# Rejecting the global T2-test means that we reject H0 along at least one
# direction, not necessarily along the coordinate direction





####### BONFERRONI INTERVALS
mu0   <- c(1,0)
alpha <- 0.01
x.mean   <- sapply(x,mean)
x.cov    <- cov(x)
x.invcov <- solve(x.cov)
n <- nrow(x) 
p <- ncol(x)  # CHANGE TO DIMENSION OF MULTIVARIETY 
k <- p        # CHANGE TO number of intervals I want to compute (set in advance)
cfr.t <- qt(1-alpha/(2*k),n-1)
Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf

# Plot of confidence and rejection regions + Bonferroni rectangle intervals
x11()
par(mfrow=c(1,1))
plot(x, asp = 1,main='Confidence and rejection regions')
ellipse(mu0, shape=x.cov/n, sqrt(cfr.fisher), col = 'blue', lty = 2, center.pch = 16)
points(x.mean[1], x.mean[2], pch = 16, col = 'red', cex=1.5)

ellipse(x.mean, shape=x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 2, center.pch = 16)
rect(Bf[1,1],Bf[2,1],Bf[1,3],Bf[2,3], border='orange', lwd=2)
legend('topleft', c('Rej. Reg.', 'Conf. Reg','Bonferroni'),col=c('blue','red','red','orange'),lty=c(2,2,1,1),lwd=2)



















####### Test (ASYMPTOTIC)
### Asymptotic test on the mean
### H0: mu == mu0 vs H1: mu != mu0
### with mu0=c(1,0)
###---------------------------------------------
# No need to verify gaussianity assumption
mu0   <- c(1,0)

x.T2A   <- n * (x.mean-mu0) %*%  x.invcov  %*% (x.mean-mu0)
cfr.chisq <- qchisq(1-alpha,p)
x.T2A < cfr.chisq # no statistical evidence to reject H0 at level alpha

# Compute the p-value
PA <- 1-pchisq(x.T2A, p)
PA


### Comparison rejection/confidence regions asymptotic or not
x11(width=14, height=7)
par(mfrow=c(1,2))
plot(x, asp = 1,main='Comparison rejection regions')
ellipse(mu0, shape=x.cov/n, sqrt(cfr.fisher), col = 'blue', lty = 1, center.pch = 4, center.cex=1.5, lwd=2)
ellipse(mu0, x.cov/n, sqrt(cfr.chisq), col = 'lightblue', lty = 1, center.pch = 4, center.cex=1.5, lwd=2)
points(mu0[1], mu0[2], pch = 4, cex = 1.5, lwd = 2, col ='lightblue')
legend('topleft', c('Exact', 'Asymptotic'),col=c('blue','lightblue'),lty=c(1),lwd=2)

plot(x, asp = 1,main='Comparison of confidence regions')
ellipse(x.mean, x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 1, center.pch = 4, center.cex=1.5, lwd=2)
ellipse(x.mean, x.cov/n, sqrt(cfr.chisq), col = 'orange', lty = 1, center.pch = 4, center.cex=1.5, lwd=2)
points(x.mean[1], x.mean[2], pch = 4, cex = 1.5, lwd = 2, col ='orange')
legend('topleft', c('Exact', 'Asymptotic'),col=c('red','orange'),lty=c(1),lwd=2)