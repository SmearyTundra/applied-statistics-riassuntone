# EXAM 15-6-2020

#############################
####### Exercise 1 ##########
#############################

pollution <- read.table('pollution.txt', header=T)

# a)
# H0: mu == c(50,50) vs H1:H0^C

# Verify the normality assumption on the data
mcshapiro.test(pollution)
# pvalue = 0.3 -> ok

n <- dim(pollution)[1]
p <- dim(pollution)[2]
mu0 <- c(50,50)
alpha <- .05

# sample mean
x.mean <- sapply(pollution,mean)
# sample cov
x.cov <- cov(pollution)
x.invcov <- solve(x.cov)

# T2 Statistics
x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
# Radius of the ellipsoid
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p) 
# Test: 
x.T2 < cfr.fisher   # FALSE 

# Compute the p-value 
P <- 1-pf(x.T2*(n-p)/((n-1)*p), p, n-p)
P 
# = 0 -> I reject H0 at level alpha


# b)
plot(pollution)
ellipse(x.mean, x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 2, lwd=2, center.cex=1)
points(mu0, col= 'blue')

# center:
x.mean

# radius
r <- sqrt(cfr.fisher)
r 
# Length of the semi-axes:
r*sqrt(eigen(x.cov/n)$values)

# c)
# The point 50,50 lies outside the confidence region for the mean, indeed we reject the null
# hypothesis of the test at point a) and we can assess that we don't have statistical evidence 
# to say that c(50,50) is the mean.

# d)
T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2

# The confidence intervals confirm my thesis: they don't containe tha value (50,50)


#############################
####### Exercise 2 ##########
#############################

stoneflakes <- read.table('stoneflakes.txt', header=T)
plot(stoneflakes)

# a)
# Euclidean distance, Ward linkage
data.e <- dist(stoneflakes, method='euclidean') # data has only the quantitative variables (other methods:
# manhattan, canberra )
data.ew <- hclust(data.e, method='ward.D2') # method: average , complete , ward.D2
# dendrogram
plot(data.ew, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(data.ew, k=3) # k=3

cluster.ew <- cutree(data.ew, k=3) 
cluster.ew

plot(stoneflakes, col = cluster.ew)

# I choose 3 clusters -> evident in the plot of my data and in dendogram!

# b)
data <- data.frame(stoneflakes, cluster.ew)
attach(data)

### Verify the assumptions:
# 1)  normality (multivariate) in each group 
Ps <- NULL
for(i in 1:3)
  Ps <- c(Ps, mcshapiro.test(data[which(cluster.ew==i),1:2])$p) 
Ps
# 0.9372 0.8404 0.4064 -> ok

# 2) same covariance structure (= same covariance matrix Sigma)
S  <-  cov(data[1:2])
S1 <-  cov(data[which(cluster.ew==1),1:2])
S2 <-  cov(data[which(cluster.ew==2),1:2])
S3 <-  cov(data[which(cluster.ew==3),1:2])


par(mfrow=c(1,3))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S3, col=heat.colors(100),main='Cov. S3', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
# ok!
dev.off()

# construct the MANOVA model:
# Model: X.ij = mu + tau.i + eps.ij; eps.ij~N_p(0,Sigma), X.ij, mu, tau.i in R^2
# Test:
# H0: tau.1 = tau.2 = tau.3  = (0,0,0)' vs H1: (H0)^c

fit <- manova(as.matrix(data[,1:2]) ~ cluster.ew)
summary.manova(fit,test="Wilks")

# I can reject the null hypothesis -> the membership to a cluster has an effect on the mean 
# features of the stone flakes! for which variable?
summary.aov(fit)

# For both variables!

# c)
alpha <- 0.1
p <- 2
g <- 3
n <- dim(stoneflakes)[1]
k <- p*g*(g-1)/2
qT <- qt(1-alpha/(2*k), n-g)

n1 <- length(which(cluster.ew==1))
n2 <- length(which(cluster.ew==2))
n3 <- length(which(cluster.ew==3))

W <- summary.manova(fit)$SS$Residuals
m  <- sapply(data[,1:2],mean)         # estimates mu
m1 <- sapply(data[which(cluster.ew==1),1:2],mean)    # estimates mu.1=mu+tau.1
m2 <- sapply(data[which(cluster.ew==2),1:2],mean)    # estimates mu.2=mu+tau.2
m3 <- sapply(data[which(cluster.ew==3),1:2],mean)    # estimates mu.3=mu+tau.3

inf12 <- m1-m2 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
sup12 <- m1-m2 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
inf13 <- m1-m3 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n3) )
sup13 <- m1-m3 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n3) )
inf23 <- m2-m3 - qT * sqrt( diag(W)/(n-g) * (1/n2+1/n3) )
sup23 <- m2-m3 + qT * sqrt( diag(W)/(n-g) * (1/n2+1/n3) )

CI <- list(cluster1.2=cbind(inf12, sup12), cluster1.3=cbind(inf13, sup13), cluster2.3=cbind(inf23, sup23))
CI

# The difference in cluster 1 and 2 in variable Length is not significant (as we can see in the
# plot for the group red and black). Otherwise, there is more difference between this two and
# the remaining one for the variable length than the variable width.
# The other intervals don't containe the value zero, so there is statistical difference in this
# groups, as we have assessed in tha MANOVA model.


#############################
####### Exercise 3 ##########
#############################

airfoil <- read.table('airfoil.txt', header=T)
attach(airfoil)

# a)
# dummy:
H <- ifelse(velocity=='H',1,0)

# model:
model <- lm(sound ~ H + frequency + frequency:H)
summary(model)

# coefficients
coeffs <- data.frame(beta0.L = model$coefficients[1],
                     beta0.H = model$coefficients[1] + model$coefficients[2],
                     beta1.L = model$coefficients[3],
                     beta1.H = model$coefficients[3] + model$coefficients[4],
                     sigma = sqrt(sum(residuals(model)^2)/model$df))
coeffs

# Verify the gaussianity assumption on the residuals
shapiro.test(model$residuals)
# pvalue = 0.52 -> ok

# Etheroschedasticity:
plot(model) # No specific patterns -> ok

# The assumptions are verified.

# b)
# 1.
model <- lm(sound ~ H + frequency + frequency:H)
linearHypothesis(model, rbind(c(0,0,1,0),
                              c(0,0,0,1)), c(0,0))
# At least one is significant. 

# 2.
linearHypothesis(model, rbind(c(1,0,0,0),
                              c(0,1,0,0)), c(0,0))
# There is statistical dependence in at least one, and looking at the summary I can assess
# that there is statistical dependence with both.

# 3.
# Looking at one at the time tests in the summary of the model,
# I can assess that the frequency with high velocity is not significant.

# c)
model.red <- lm(sound ~ H + frequency)
summary(model.red)

coeffs.new <- data.frame(beta0.L = model.red$coefficients[1],
                         beta0.H = model.red$coefficients[1] + model.red$coefficients[2],
                         beta1.L = model.red$coefficients[3],
                         sigma = sqrt(sum(residuals(model.red)^2)/model.red$df))
coeffs.new

# d)
x.new <- data.frame(H = 1, frequency = 15000)
predict(model.red,x.new, interval = 'confidence')
#      fit      lwr      upr
# 417.6368 393.8384 441.4352

