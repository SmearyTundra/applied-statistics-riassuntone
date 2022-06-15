#########
### Tests and confidence regions for the mean of a multivariate Gaussian
########

library(car)
library(mvtnorm)

stiff <- read.table('', header = T)
head(stiff)
dim(stiff)

n <- dim(stiff)[1]
p <- dim(stiff)[2]

plot(stiff,pch=19)

dev.off()

### test of Gaussianity
mcshapiro.test(stiff)

###----------------------------
### Test for the mean of level 5%
###--------------------------------
# 1) Formulate the test:
#    H0: mu == mu0 vs H1: mu != mu0
#    with mu0=c(1850, 1750, 1500, 1700)

mu0   <- c(1850, 1750, 1500, 1700)
alpha <- 0.05

# 2) Compute the test statistics
x.mean   <- colMeans(stiff)
x.cov    <- cov(stiff)
x.invcov <- solve(x.cov)

x.T2 <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 

# 3a) Verify if the test statistics belongs to the rejection region
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)                     #If only one variable  qchisq(0.9,p)
x.T2 < cfr.fisher # we accept H0 if TRUE (in this case at lever 5%)

# 3b) Compute the p-value
P <- 1-pf(x.T2*(n-p)/((n-1)*p), p, n-p)
P

###------------------------------------------------------------
### Confidence region for the mean of level 95%
###------------------------------------------------------------
#CR for the mean (ellipsoidal region) 
#     { m \in R^4(?) t.c. n * (x.mean-m)' %*% (x.cov)^-1 %*% (x.mean-m) < cfr.fisher }

### test of Gaussianity
mcshapiro.test(stiff)

# Center:
x.mean

# Directions of the principal axes:
eigen(x.cov/n)$vectors

# Radius
r <- sqrt(cfr.fisher)

# Length of the semi-axes of the ellipse:
r*sqrt(eigen(x.cov/n)$values) 

#ellipse
plot(stiff,pch=19)
ellipse(x.mean, x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 2, lwd=2, center.cex=1)
points(x.mean[1], x.mean[2], pch = 16, col ='red', cex = 1.5)

# Region of rejection (centered in mu0)
plot(stiff,pch=19)
ellipse(mu0, shape=x.cov/n, sqrt(cfr.fisher), col = 'blue', lty = 2, center.pch = 16)
points(x.mean[1], x.mean[2], pch = 16, col ='red', cex = 1.5)

# Remark: the radius and the shape of the ellipse are the same, but the center changes:
# - Rejection region: the center is the mean mu0 under H0 (blue ellipse)
# - Confidence region: the center is the sample mean (red ellipse)

#####
# Which relation between the two ellipses?
# - If the rejection region does NOT contain the sample mean (i.e., we
#   are in the acceptance region), then we cannot reject H0 
#   (i.e., if the sample mean falls within the ellipse we accept H0)
# - If the mean under H0 (mu0) is contained in the confidence region
#   of level 1-alpha, then we do not reject H0 at level alpha
# => the confidence region of level 1-alpha contains all the mu0
#    that we would accept at level alpha

###------------------------------------------------------------
### Confidence region for 2 independent population with n1!=n2
###------------------------------------------------------------
#CR for the mean (ellipsoidal region) 
#     { m \in R^4(?) t.c. n * (x.mean-m)' %*% (x.cov)^-1 %*% (x.mean-m) < cfr.fisher }

n1 <- length(c1) 
n2 <- length(c2)
p <- dim(whales)[2]

m1 <- colMeans(whales[c1,])
m2 <- colMeans(whales[c2,])

mcshapiro.test(whales)
shapiro.test(whales[c1,2])$p.value
shapiro.test(whales[c2,2])$p.value

S1 <- cov(whales[c1,])
S2 <- cov(whales[c2,])
Sp <- (S1 * (n1 - 1) + S2 * (n2 - 1)) / (n1 + n2 -2)
invSP <- solve(Sp)

# Center:
m1 -m2

# Directions of the principal axes:
eigen(Sp)$vectors  

# Radius
cfr.fisher <- p*(n1 + n2 - 2)/(n1 + n2 - 1 - p)*qf(1-alpha, p, n1 + n2 - 1 - p)
r <- sqrt(cfr.fisher)  

# Length of the semi-axes of the ellipse:
r*sqrt(eigen(Sp*(1/n1+1/n2))$values)  

plot(m1[1]-m2[1], m1[2]-m2[2], asp = 1, xlim = c(8,13))
ellipse(m1 - m2, shape=Sp*(1/n1+1/n2), sqrt(cfr.fisher), col = 'blue', lty = 2, center.pch = 20)  


###----------------------------------------------------
### Simultaneous T2 intervals on the components of the mean with global level 95%
###----------------------------------------------------

T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2

###----------------------------------------------------
### Bonferroni intervals on the components of the mean with global level 95%
###----------------------------------------------------
k <- p
cfr.t <- qt(1 - alpha/(k*2), n-1)

Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf

###
# Let's do a plot
matplot(1:4,1:4,pch='',ylim=range(stiff),xlab='Variables',ylab='Confidence intervals along a component',main='Confidence intervals')

for(i in 1:4) segments(i,T2[i,1],i,T2[i,3],lwd=2,col='grey35', lty=3)
points(1:4, T2[,1], pch='-', col='grey35')
points(1:4, T2[,3], pch='-', col='grey35')

for(i in 1:4) segments(i,Bf[i,1],i,Bf[i,3],lwd=2,col=i)
points(1:4, Bf[,2], pch=16, col=1:4)
points(1:4, Bf[,1], pch='-', col=1:4)
points(1:4, Bf[,3], pch='-', col=1:4)

# Is mu0 inside the Bonferroni confidence region?
# we add it to the plot
points(1:4, mu0, lwd=3, col='orange')





# Estimate an elliptical region A that contains 99% of the pines
# Ellipse(mean,S,r)
# And not S/n as you do to find the elliptical confidence region for the mean