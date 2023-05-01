remove(list = setdiff(ls(), lsf.str()))
library(car)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#
#         CONFIDENCE INTERVALS FOR GAUSSIAN MEAN AND HYPOTHESIS TESTING
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# For some dataframe x
# Setting up
# Dimensions
n <- dim(x)[1]
p <- dim(x)[2]
alpha <- 0.05 # At level 5%

# Sample mean, covariance 
x.mean   <- sapply(x,mean)
x.cov    <- cov(x)
x.invcov <- solve(x.cov)
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)

# T2 simultaneous intervals

T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2

# Bonferroni confidence intervals

k <- p # number of intervals I want to compute (set in advance)
cfr.t <- qt(1-alpha/(2*k),n-1)
Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf

# If you want to put everything in a loop for bonferroni

alpha <- 0.1 # At level 5%
k <- 10
#cluster
# For means within individual clusters

for(i in 1:5) {
  # Sample mean, covariance 
  x.mean   <- sapply(sequoia[get(paste('i',i, sep='')),1:2],mean)
  x.cov    <- cov(sequoia[get(paste('i',i, sep='')),1:2])
  x.invcov <- solve(x.cov)
  cfr.t <- qt(1-alpha/(2*k),n-1)
  print(paste("Mean of cluster",i))
  ICMean <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
                  center = x.mean, 
                  sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
  # We perform a chi-square test on the variance
  ICVar <- cbind(inf=diag(x.cov)*(n-1) / qchisq(1 - alpha/(2*k), n-1),
                 center=diag(x.cov),
                 sup=diag(x.cov)*(n-1) / qchisq(alpha/(2*k), n-1))
  print(ICMean)
  print(paste("Variance of cluster",i))
  print(ICVar)
}

# For some mu0 you can perform the test H0: mu == mu0

x.T2 <- n * t(x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
x.T2 < cfr.fisher

# Compute the p-value
P <- 1-pf(x.T2*(n-p)/((n-1)*p), p, n-p)
P
                
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                CONFIDENCE INTERVALS FOR DIFFERENCE 
#                GAUSSIAN MEAN AND HYPOTHESIS TESTING
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Assumption EQUAL VARIANCES!!!
# For some dataframes t1 and t2

n1 <- dim(t1)[1] #
n2 <- dim(t2)[1] #
p  <- dim(t1)[2] # 

# Compute sample means and sample covariances
t1.mean <- sapply(t1,mean)
t2.mean <- sapply(t2,mean)
t1.cov  <-  cov(t1)
t2.cov  <-  cov(t2)

# Compute weighted covariance Spooled
Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2)
# we compare the matrices
list(S1=t1.cov, S2=t2.cov, Spooled=Sp)

# Simultaneous T2 intervals (assuming 2 dimensions!!!)
IC.T2.X1 <- c(t1.mean[1]-t2.mean[1]-sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)), 
              t1.mean[1]-t2.mean[1]+sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)) )
IC.T2.X2 <- c(t1.mean[2]-t2.mean[2]-sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)), 
              t1.mean[2]-t2.mean[2]+sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)) )
IC.T2 <- rbind(IC.T2.X1, IC.T2.X2)
dimnames(IC.T2)[[2]] <- c('inf','sup')                        
IC.T2

# Bonferroni intervals

k <- p # number of intervals I want to compute (set in advance)
cfr.t <- qt(1-alpha/(2*k),n1+n2-2)
Bf1 <- cbind(inf = (t1.mean[1]-t2.mean[1]) - cfr.t*sqrt(Sp[1,1]*(1/n1+1/n2)),
            sup = x(t1.mean[1]-t2.mean[1]) + cfr.t*sqrt(Sp[1,1]*(1/n1+1/n2)))
Bf2 <- cbind(inf = (t1.mean[2]-t2.mean[2]) - cfr.t*sqrt(Sp[2,2]*(1/n1+1/n2)),
             sup = x(t1.mean[2]-t2.mean[2]) + cfr.t*sqrt(Sp[2,2]*(1/n1+1/n2)))
Bf <- rbind(Bf1, Bf2)
dimnames(Bf)[[2]] <- c('inf','sup')    
Bf

# HYPOTHESIS TESTING
# For some mu0

# Test H0: mu1 == mu2  vs  H1: mu1 != mu2
# i.e.,
# Test H0: mu1-mu2 == c(0,0)  vs  H1: mu1-mu2 != c(0,0)

alpha   <- .01
delta.0 <- c(0,0)
Spinv   <- solve(Sp)

T2 <- n1*n2/(n1+n2) * (t1.mean-t2.mean-delta.0) %*% Spinv %*% (t1.mean-t2.mean-delta.0)

cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
T2 < cfr.fisher # TRUE: no statistical evidence to reject H0 at level 1%

P <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
P  


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#               ELLIPTICAL REGION FOR GAUSSIAN MEAN, SEMIAXES
#               CENTER, LENGTH OF AXES, ETC
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

n <- dim(D)[1]  # 11
p <- dim(D)[2]  #  2
D.mean   <- sapply(D,mean)
D.cov    <- cov(D)
D.invcov <- solve(D.cov)
alpha   <- .05

# Distributed as Fisher with p, n-p degrees of freedom
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
cfr.fisher
# Spectral decomposition of convariance matrix
decomp <- eigen(D.cov)
decomp
# Direction of axes
decomp$vectors
# Center
D.mean
# Radius of ellipse
r <- sqrt(cfr.fisher)
r
#Length of semi-axes
lengthSemiAxes <- r*sqrt(decomp$values)
lengthSemiAxes


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                                 PAIRED COMPARISONS
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Compute differences (for mean mu0 = c(0,0))
# WATCH OUT FOR THE DIFFERENCE YOU WANT!
D <- data.frame(diff1=dataset[,1]-dataset[,3], diff2=dataset[,2]-dataset[,4]) 
D

x11()
plot(D, asp=1, pch=19, main='Dataset of Differences')
abline(h=0, v=0, col='grey35')
points(0,0, pch=19, col='grey35')

# Multivariate shapiro on D
mcshapiro.test(D)

# Computation of Hotelling's T2 and quantile of Fisher distribution

n <- dim(D)[1]  # 11
p <- dim(D)[2]  #  2

D.mean   <- sapply(D,mean)
D.cov    <- cov(D)
D.invcov <- solve(D.cov)

# Significance level

alpha   <- .05
delta.0 <- c(0,0)

# Hotelling's T2
D.T2 <- n * (D.mean-delta.0) %*% D.invcov %*% (D.mean-delta.0)
D.T2

# Fisher quantile
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
cfr.fisher

# Test result
D.T2 < cfr.fisher # FALSE: we reject H0 at level 5%

# P-value computation
P <- 1-pf(D.T2*(n-p)/(p*(n-1)), p, n-p)
P

# Ellipsoidal confidence region with confidence level (1-alpha)100%
x11()
plot(D, asp=1, pch=1, main='Dataset of the Differences',ylim=c(-15,60))

# Ellipse centered around our sample mean
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt(cfr.fisher), lwd=2)

# Ellipse centered around mu0
ellipse(center=delta.0, shape=D.cov/n, radius=sqrt(cfr.fisher), lwd=2, color="blue")
points(delta.0[1], delta.0[2], pch=16, col='red', cex=1.5)
abline(h=delta.0[1], v=delta.0[2], col='black')
abline(h=0,v=0)
dev.off()

# Simultatenous T2 intervals
IC.T2.DBOD <- c( D.mean[1]-sqrt(cfr.fisher*D.cov[1,1]/n) , D.mean[1], D.mean[1]+sqrt(cfr.fisher*D.cov[1,1]/n) )
IC.T2.DSS  <- c( D.mean[2]-sqrt(cfr.fisher*D.cov[2,2]/n) , D.mean[2], D.mean[2]+sqrt(cfr.fisher*D.cov[2,2]/n) )

T2 <- rbind(IC.T2.DBOD, IC.T2.DSS)
dimnames(T2)[[2]] <- c('inf','center','sup')
T2

# Enough to reject along one direction, find such direction (worst)

D.T2
#    - the distribution of the maximum is known
#    - the direction along which the maximum is realized is known
worst <- D.invcov %*% (D.mean-delta.0)
worst <- worst/sqrt(sum(worst^2))
worst
# Angle with the x-axis:
theta.worst <- atan(worst[2]/worst[1])+pi
theta.worst

# Confidence interval along the worst direction:
IC.worst  <- c( D.mean %*% worst - sqrt(cfr.fisher*(t(worst)%*%D.cov%*%worst)/n),
                D.mean %*% worst,
                D.mean %*% worst + sqrt(cfr.fisher*(t(worst)%*%D.cov%*%worst)/n) )
IC.worst
delta.0%*%worst
(IC.worst[1] < delta.0%*%worst) & (delta.0%*%worst < IC.worst[2])   
# Reject H0: a'mu == a'delta.0 in direction a=worst

# Extremes of IC.worst in the coordinate system (x,y):
x.min <- IC.worst[1]*worst
x.max <- IC.worst[3]*worst
m1.ort <- -worst[1]/worst[2]
q.min.ort <- x.min[2] - m1.ort*x.min[1]
q.max.ort <- x.max[2] - m1.ort*x.max[1]
abline(q.min.ort, m1.ort, col='forestgreen', lty=2,lwd=1)
abline(q.max.ort, m1.ort, col='forestgreen', lty=2,lwd=1)

m1=worst[2]/worst[1] # worst direction
abline(0, m1, col='grey35')
segments(x.min[1],x.min[2],x.max[1],x.max[2],lty=1,lwd=2, col='forestgreen')

dev.off()

# Bonferroni intervals for the mean
k <- p  # number of parameters or as many as you want (for k directions)
cfr.t <- qt(1-alpha/(2*k),n-1)

IC.BF.DBOD <- c( D.mean[1]-cfr.t*sqrt(D.cov[1,1]/n) , D.mean[1], D.mean[1]+cfr.t*sqrt(D.cov[1,1]/n) )
IC.BF.DSS  <- c( D.mean[2]-cfr.t*sqrt(D.cov[2,2]/n) , D.mean[2], D.mean[2]+cfr.t*sqrt(D.cov[2,2]/n) )

Bf <- rbind(IC.BF.DBOD, IC.BF.DSS)
dimnames(Bf)[[2]] <- c('inf','center','sup')
Bf

x11()
plot(D, asp=1, pch=1, main='Dataset of the Differences',ylim=c(-15,60))
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt(cfr.fisher), lwd=2, col='grey', center.cex=1.25)

abline(h=0, v=0, col='grey', lty=1, lwd=2)
points(delta.0[1], delta.0[2], pch=16, col='grey35', cex=1.25)

abline(v = T2[1,1], col='red', lwd=1, lty=2)
abline(v = T2[1,3], col='red', lwd=1, lty=2)
abline(h = T2[2,1], col='red', lwd=1, lty=2)
abline(h = T2[2,3], col='red', lwd=1, lty=2)
segments(IC.T2.DBOD[1],0,IC.T2.DBOD[3],0,lty=1,lwd=2,col='red')
segments(0,IC.T2.DSS[1],0,IC.T2.DSS[3],lty=1,lwd=2,col='red')

abline(v = Bf[1,1], col='blue', lwd=1, lty=2)
abline(v = Bf[1,3], col='blue', lwd=1, lty=2)
abline(h = Bf[2,1], col='blue', lwd=1, lty=2)
abline(h = Bf[2,3], col='blue', lwd=1, lty=2)
segments(IC.BF.DBOD[1],0,IC.BF.DBOD[3],0,lty=1,lwd=2,col='blue')
segments(0,IC.BF.DSS[1],0,IC.BF.DSS[3],lty=1,lwd=2,col='blue')
dev.off()


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                        TESTS FOR REPEATED MEASURES
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Tests for repeated measures (same thing measured multiple times under different circumstances)

# IF the x axis can be seen as time then it makes sense to make a matplot
# For some dataset 

matplot(t(dataset), type='l')
n <- dim(dataset)[1]
q <- dim(dataset)[2]

M <- sapply(dataset,mean)
M
S <- cov(dataset)
S

# Contrast matrix, build it depending on what comparisons you want to later. For 
# hypothesis testing only any contrast matrix will do
C <- matrix(c(-1, 1, 0, 0,
              -1, 0, 1, 0,
              -1, 0, 0, 1), 3, 4, byrow=T)
C

# Computation of our statistics
# Test: H0: C%*%mu == 0 vs H1: C%*%mu != 0
alpha   <- .05
delta.0 <- c(0,0,0)

Md <- C %*% M 
Sd <- C %*% S %*% t(C)
Sdinv <- solve(Sd)

# Hotelling's T2

T2 <- n * t( Md - delta.0 ) %*% Sdinv %*% ( Md - delta.0 )

cfr.fisher <- ((q-1)*(n-1)/(n-(q-1)))*qf(1-alpha,(q-1),n-(q-1)) 

# Reject or don't reject H0

T2 < cfr.fisher
T2

cfr.fisher

# T2 is much higher than cfr.fisher => the p-value will be very small

P <- 1-pf(T2*(n-(q-1))/((q-1)*(n-1)),(q-1),n-(q-1))
P

# T2 simultaneous intervals for the differences

IC.T2 <- cbind( Md - sqrt(cfr.fisher*diag(Sd)/n) , Md, Md + sqrt(cfr.fisher*diag(Sd)/n) )
IC.T2

# Bonferroni intervals for the differences

k     <- q - 1   # number of increments (i.e., dim(C)[1]) (levels whatever)
cfr.t <- qt(1-alpha/(2*k),n-1)
IC.BF <- cbind( Md - cfr.t*sqrt(diag(Sd)/n) , Md, Md + cfr.t*sqrt(diag(Sd)/n) )
IC.BF

# TEST FOR mu0 OTHER THAN ZERO

### what if we want to verify the following hypothesis:
### "the drug decreases the pressure of two units with respect to
### the baseline at both 8 and 16 hours, and its effect vanishes in 24 hours
### from the drug administration"

C <- matrix(c(-1, 1, 0, 0,
              -1, 0, 1, 0,
              -1, 0, 0, 1), 3, 4, byrow=T)
delta.0 <- c(-2,-2,0)

# or

C <- matrix(c(-1, 1, 0, 0,
              0, -1, 1, 0,
              0, 0, -1, 1), 3, 4, byrow=T)
delta.0 <- c(-2,0,2)

# To visualize intervals
x11()
matplot(t(matrix(1:3,3,3)),t(IC.BF), type='b',pch='',xlim=c(0,4),xlab='',ylab='', main='Confidence intervals')
segments(matrix(1:3,3,1),IC.BF[,1],matrix(1:3,3,1),IC.BF[,3], col='orange', lwd=2)
points(1:3, IC.BF[,2], col='orange', pch=16)
points(1:3+.05, delta.0, col='black', pch=16)
segments(matrix(1:3+.1,3,1),IC.T2[,1],matrix(1:3+.1,3,1),IC.T2[,3], col='blue', lwd=2)
points(1:3+.1,IC.T2[,2], col='blue', pch=16)
legend('topright', c('Bonf. IC', 'Sim-T2 IC'), col=c('orange', 'blue'), lty=1, lwd=2)





# TEST FOR THE MEAN OF TWO INDEPENDENT GAUSSIAN POPULATIONS
# Covariance matrix estimate with Spooled

# we build the data
t1 <- matrix(c(3,3,1,6,2,3),2)
t1 <- data.frame(t(t1))
t2 <- matrix(c(2,3,5,1,3,1,2,3),2)
t2 <- data.frame(t(t2))

t1
t2

n1 <- dim(t1)[1] # n1=3
n2 <- dim(t2)[1] # n2=4
p  <- dim(t1)[2] # p=2


# we compute the sample mean, covariance matrices and the matrix Spooled

t1.mean <- sapply(t1,mean)
t2.mean <- sapply(t2,mean)
t1.cov  <-  cov(t1)
t2.cov  <-  cov(t2)
Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2)
# we compare the matrices
list(S1=t1.cov, S2=t2.cov, Spooled=Sp)

# Test H0: mu1 == mu2  vs  H1: mu1 != mu2
# i.e.,
# Test H0: mu1-mu2 == c(0,0)  vs  H1: mu1-mu2 != c(0,0)

alpha   <- .01
delta.0 <- c(0,0)
Spinv   <- solve(Sp)

T2 <- n1*n2/(n1+n2) * (t1.mean-t2.mean-delta.0) %*% Spinv %*% (t1.mean-t2.mean-delta.0)

cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
T2 < cfr.fisher # TRUE: no statistical evidence to reject H0 at level 1%

P <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
P  
# P-value high (we don't reject at 1%,5%,10%)

# Simultaneous T2 intervals
IC.T2.X1 <- c(t1.mean[1]-t2.mean[1]-sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)), t1.mean[1]-t2.mean[1]+sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)) )
IC.T2.X2 <- c(t1.mean[2]-t2.mean[2]-sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)), t1.mean[2]-t2.mean[2]+sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)) )
IC.T2 <- rbind(IC.T2.X1, IC.T2.X2)
dimnames(IC.T2)[[2]] <- c('inf','sup')                        
IC.T2

graphics.off()

# FOR BONFERRONI
x.mean   <- colMeans(data)
x.cov    <- cov(data)
x.invcov <- solve(x.cov)
k <- 2 # number of intervals I want to compute (set in advance)
cfr.t <- qt(1-alpha/(2*k),n-1)
Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                             FALSE DISCOVERY RATE
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# You'll have to correct the univariate p-values of multiple tests.
# R function to correct p-values: p.adjust

allergy <- read.table('hatingalmonds.txt')
head(allergy)
dim(allergy)

noallergy <- read.table('lovingalmonds.txt')
head(noallergy)
dim(noallergy)

n1 <- dim(allergy)[1]
n2 <- dim(noallergy)[1]
p <- dim(noallergy)[2]

x.mean1 <- sapply(allergy, mean)
x.mean2 <- sapply(noallergy, mean)

p.hat <- (x.mean1*n1+x.mean2*n2)/(n1+n2)
x.var <- (p.hat*(1-p.hat))

# Test: H0.i: mu.i1 == mu.i2  vs H1.i: mu.i1 != mu.i2

z.i <- (x.mean1-x.mean2)/sqrt(x.var*(1/n1+1/n2))
p.i <- ifelse(z.i<0, 2*pnorm(z.i),2*(1-pnorm(z.i)))

which(p.i<.01)

# Bonferoni test
k <- 520

which(p.i*k<.01)  

# or
p.Bf <- p.adjust(p.i, method='bonferroni')

which(p.Bf<.01)  

# Benjamini-Hochberg (control the false discovery rate)  
p.BH <- p.adjust(p.i, method='BH')

which(p.BH<.01)


x11(width=21, height=7)
par(mfrow=c(1,3))
plot(p.i, main='Univariate')
abline(h=.01, lwd=2, col='red')

plot(p.Bf, main='Corrected - Bonferroni')
abline(h=.01, lwd=2, col='red')

plot(p.BH, main='Corrected - BH')
abline(h=.01, lwd=2, col='red')


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                        UNIVARIATE T-TEST GAUSSIAN MEAN
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# One-sided
# Consider the fourth column, the mean of that column must be greater than 0.2 times the mean of the first one
# Comparison of the means of two gaussian populations
# t-test from Stat 101
# H0: mu1 >= mu0*0.2 
# H0: - mu1 + mu0*0.2 <= 0
# H0: a'mu <= 0 vs H1: a'mu < 0 where a = c(0.2,-1)

a <- as.matrix(c(0.2,-1), 2,1)
delta.0 <- 0
accPurchase <- as.matrix(cbind(D[,1],D[,4]))
# Must perform a t-test
t.stat <- (mean(accPurchase %*% a) - delta.0) / sqrt( var(accPurchase %*% a) / n ) 
t.stat

P <- 1-pt(t.stat, n-1)
P

# VANILLA TEST FOR MULTIVARIATE MEAN OF GAUSSIAN

n <- dim(x)[1]
p <- dim(x)[2]

x.mean   <- sapply(x,mean)
x.cov    <- cov(x)
x.invcov <- solve(x.cov)

alpha <- 0.01
mu0 <- c(1,0)

# T2 SIMULTANTEOUS INTERVALS!
# T2 Statistics
x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
# Radius of the ellipsoid
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)

T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2

x11()
matplot(1:4,1:4,pch='',ylim=range(stiff),xlab='Variables',ylab='T2 for a component', 
        main='Simultaneous T2 conf. int. for the components')
for(i in 1:4)segments(i,T2[i,1],i,T2[i,3],lwd=3,col=i)
points(1:4, T2[,2], pch=16, col=1:4)

# BONFERRONI INTERVALS

k <- p
cfr.t <- qt(1 - alpha/(k*2), n-1)

Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf

# PLOT OF T2 INTERVALS
x11()
matplot(1:4,1:4,pch='',ylim=range(dataset),xlab='Variables',ylab='T2 for a component', 
        main='Simultaneous T2 conf. int. for the components')
for(i in 1:4)segments(i,T2[i,1],i,T2[i,3],lwd=3,col=i)
points(1:4, T2[,2], pch=16, col=1:4)


# FOR CLASSIC MEAN CONFIDENCE REGION


### Simultaneous T2 intervals on the components of the mean
### with global level 95%
###----------------------------------------------------

T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2

x11()
matplot(1:4,1:4,pch='',ylim=range(stiff),xlab='Variables',ylab='T2 for a component', 
        main='Simultaneous T2 conf. int. for the components')
for(i in 1:4)segments(i,T2[i,1],i,T2[i,3],lwd=3,col=i)
points(1:4, T2[,2], pch=16, col=1:4)

# Is mu0 inside the rectangular region?
# We add it to the plot
points(1:4, mu0, lwd=3, col='orange')

# Yes, it is, because it is inside all the T2-intervals,

### Bonferroni intervals on the components of the mean
### with global level 95%
###----------------------------------------------------
k <- p
cfr.t <- qt(1 - alpha/(k*2), n-1)

Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf

# Let's do a plot
x11()
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

# Yes, it is, because it belongs to all the intervals along the components

graphics.off()


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#
#               Stat 101 - Mean and variance univariate case, plus Bonferroni
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Consider a dataset d
x <- c(1,2,3,4,5,1,5,5,6,6,2,2,2,22,2,2,23,3,2,2,2,2,2)

x.mean <- mean(x)
x.var <- var(x)
n <- length(x)

alpha <- 0.05
k <- 1
# t-student quantile (unknown variance)

cfr.t <- qt(1 - alpha/(2*k), n-1)

# Our confidence interval for the mean

Bfmean <- cbind(inf = x.mean - cfr.t*sqrt(x.var/n),
            center = x.mean, 
            sup = x.mean + cfr.t*x.var*sqrt(x.var/n))
Bfmean


# Confidence interval for the variance (chi square distribution)

Bfvar <- cbind(inf=x.var*(n-1) / qchisq(1 - alpha/(2*k), n-1),
               center=x.var,
               sup=x.var*(n-1) / qchisq(alpha/(2*k), n-1))
Bfvar

# Want Bonferroni intervals? Divide alpha by k and that's it!