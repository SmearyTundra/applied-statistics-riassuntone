### TOPICS:
### Test for the mean of paired multivariate Gaussian observations
### Test for repeated measures
### Test for two independent Gaussian populations
library(car)


#### TO ADD 
## FDR (FILE GRUPPO) 





############# TEST MEAN PAIRED MULTIVARIATE GAUSSIAN ##############
# Compute differences (for mean mu0 = c(0,0))
# WATCH OUT FOR THE DIFFERENCE YOU WANT!
d <- data.frame(diff1=data[,1]-data[,3], diff2=data[,2]-data[,4]) 

x11()
plot(d, asp=1, pch=19, main='data of differences')
abline(h=0, v=0, col='grey35')
points(0,0, pch=19, col='grey35')

# Multivariate shapiro on d
mcshapiro.test(d)


#### TEST
# Computation of Hotelling's T2 and quantile of Fisher distribution
n <- dim(d)[1]  # 11
p <- dim(d)[2]  #  2
d.mean   <- sapply(d,mean)
d.cov    <- cov(d)
d.invcov <- solve(d.cov)

# Significance level
alpha   <- .05
delta.0 <- c(0,0) # equivalent mu0

# Hotelling's T2
d.T2 <- n * (d.mean-delta.0) %*% d.invcov %*% (d.mean-delta.0)
d.T2
# Fisher quantile
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
cfr.fisher

# Test result
d.T2 < cfr.fisher # FALSE: we reject H0 at level 5%

# P-value computation
P <- 1-pf(d.T2*(n-p)/(p*(n-1)), p, n-p)
P


#### CONFIDENCE/REJECTION REGION
# Ellipsoidal confidence region with confidence level (1-alpha)100%
x11()
plot(d, asp=1, pch=1, main='data of the differences',ylim=c(-15,60))

# Ellipse centered around our sample mean
ellipse(center=d.mean, shape=d.cov/n, radius=sqrt(cfr.fisher), lwd=2)

# Ellipse centered around mu0 -> REJECTION
ellipse(center=delta.0, shape=d.cov/n, radius=sqrt(cfr.fisher), lwd=2, color="red")
points(delta.0[1], delta.0[2], pch=16, col=1, cex=1.5)
abline(h=delta.0[1], v=delta.0[2], col='black')
abline(h=0,v=0)
dev.off()




#### SIMULTANEOUS T2 INTERVALS FOR THE MEAN
IC.T2.1 <- c( d.mean[1]-sqrt(cfr.fisher*d.cov[1,1]/n) , d.mean[1], d.mean[1]+sqrt(cfr.fisher*d.cov[1,1]/n) )
IC.T2.2  <- c( d.mean[2]-sqrt(cfr.fisher*d.cov[2,2]/n) , d.mean[2], d.mean[2]+sqrt(cfr.fisher*d.cov[2,2]/n) )
T2 <- rbind(IC.T2.1, IC.T2.2)
dimnames(T2)[[2]] <- c('inf','center','sup')
T2                
# plot of T2 intervalse and confidence region
x11()
plot(d, asp=1, pch=1, main='dataset of the differences',ylim=c(-15,60))
ellipse(center=d.mean, shape=d.cov/n, radius=sqrt(cfr.fisher), lwd=2, col='grey')
abline(v = T2[1,1], col='red', lwd=1, lty=2)
abline(v = T2[1,3], col='red', lwd=1, lty=2)
abline(h = T2[2,1], col='red', lwd=1, lty=2)
abline(h = T2[2,3], col='red', lwd=1, lty=2)
points(delta.0[1], delta.0[2], pch=16, col='grey35', cex=1.5)
abline(h=delta.0[1], v=delta.0[2], col='grey35')


#### WORST DIRECTION ANALYSIS
# Enough to reject along one direction, find such direction (worst)
#    - the distribution of the maximum is known
#    - the direction along which the maximum is realized is known
worst <- d.invcov %*% (d.mean-delta.0)
worst <- worst/sqrt(sum(worst^2))
worst
# Angle with the x-axis:
theta.worst <- atan(worst[2]/worst[1])+pi
theta.worst

# Confidence interval along the worst direction:
IC.worst  <- c( d.mean %*% worst - sqrt(cfr.fisher*(t(worst)%*%d.cov%*%worst)/n),
                d.mean %*% worst,
                d.mean %*% worst + sqrt(cfr.fisher*(t(worst)%*%d.cov%*%worst)/n) )
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





#### BONFERRONI INTERVALS FOR THE MEAN
n <- dim(d)[1]  # 11
p <- dim(d)[2]  #  2
d.mean   <- sapply(d,mean)
d.cov    <- cov(d)
d.invcov <- solve(d.cov)
k <- p  # number of parameters or as many as you want (for k directions)

# Significance level
alpha   <- .05
delta.0 <- c(0,0) # equivalent mu0
d.T2 <- n * (d.mean-delta.0) %*% d.invcov %*% (d.mean-delta.0)
cfr.t <- qt(1-alpha/(2*k),n-1)

IC.BF.1 <- c( d.mean[1]-cfr.t*sqrt(d.cov[1,1]/n) , d.mean[1], d.mean[1]+cfr.t*sqrt(d.cov[1,1]/n) )
IC.BF.2  <- c( d.mean[2]-cfr.t*sqrt(d.cov[2,2]/n) , d.mean[2], d.mean[2]+cfr.t*sqrt(d.cov[2,2]/n) )

Bf <- rbind(IC.BF.1, IC.BF.2)
dimnames(Bf)[[2]] <- c('inf','center','sup')
Bf


ICvar <- cbind(inf     = diag(d.cov)*(n-1) / qchisq(1 - alpha/(2*k), n-1),
           center  = diag(d.cov),
           sup     = diag(d.cov)*(n-1) / qchisq(alpha/(2*k), n-1))
ICvar



# Plot of confidence + T2 + bonf
x11()
plot(d, asp=1, pch=1, main='data of the differences',ylim=c(-15,60))
ellipse(center=d.mean, shape=d.cov/n, radius=sqrt(cfr.fisher), lwd=2, col='grey', center.cex=1.25)
abline(h=0, v=0, col='grey', lty=1, lwd=2)
points(delta.0[1], delta.0[2], pch=16, col='grey35', cex=1.25)

abline(v = T2[1,1], col='red', lwd=1, lty=2)
abline(v = T2[1,3], col='red', lwd=1, lty=2)
abline(h = T2[2,1], col='red', lwd=1, lty=2)
abline(h = T2[2,3], col='red', lwd=1, lty=2)

abline(v = Bf[1,1], col='blue', lwd=1, lty=2)
abline(v = Bf[1,3], col='blue', lwd=1, lty=2)
abline(h = Bf[2,1], col='blue', lwd=1, lty=2)
abline(h = Bf[2,3], col='blue', lwd=1, lty=2)
dev.off()


# 2) Bonf. for means (ALTERNATIVE)
x.mean   <- colMeans(data)
x.cov    <- cov(data)
x.invcov <- solve(x.cov)
k <- 2 # number of intervals I want to compute (set in advance)
cfr.t <- qt(1-alpha/(2*k),n-1)
Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf



















############# TEST MEAN TWO INDEPENDENT MULTIVARIATE GAUSSIAN ##############
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



#### TEST
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
P   # P-value high (we don't reject at 1%,5%,10%)



#### SIMULTANEOUS T2 INTERVALS FOR THE DIFF. oF MEANS
IC.T2.X1 <- c(t1.mean[1]-t2.mean[1]-sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)), t1.mean[1]-t2.mean[1]+sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)) )
IC.T2.X2 <- c(t1.mean[2]-t2.mean[2]-sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)), t1.mean[2]-t2.mean[2]+sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)) )
IC.T2 <- rbind(IC.T2.X1, IC.T2.X2)
dimnames(IC.T2)[[2]] <- c('inf','sup')                        
IC.T2



#### BONFERRONI INTERVALS FOR THE DIFFERENCE OF MEANS
# 1) Bonf intervals for the difference of means
alpha <- 0.1
IC <- cbind(t2.mean-t1.mean - sqrt(diag(Sp)*(1/n1+1/n2)) * qt(1 - alpha/(p*2), n1+n2-2),
            t2.mean-t1.mean,
            t2.mean-t1.mean + sqrt(diag(Sp)*(1/n1+1/n2)) * qt(1 - alpha/(p*2), n1+n2-2))
IC




#### BONUS: MANOVA TEST
S1 <- cov(girona)
S2 <- cov(terrassa)
x11(width=21)
par(mfrow=c(1,2))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
##le varianze possono essere considerate abbastanza simili
### posso procedere con le analisi

#USANDO MANOVA
smean_girona <- sapply(girona, mean)
smean_terrassa <- sapply(terrassa, mean)

data <- rbind(girona, terrassa)
fact <- as.factor(c(rep("girona", 35), rep("terrassa", 35)))

#we perform a manova test
man <- manova(as.matrix(data) ~ fact)
summary.aov(man)
summary.manova(man)
#we have statistical evidence to reject the null hypothesis H0: mu_girona = mu_terrassa
mcshapiro.test(man$residuals)

























############# REPEATED MEASURES ##############

d <- read.table ('d.txt', col.names=c('h.0','h.8','h.16','h.24'))

mcshapiro.test(d)
x11()
matplot(t(d), type='l')

### Perform a test at level 5% to prove that there is a difference
### between the means of repeated measures
n <- dim(d)[1]
q <- dim(d)[2]

M <- sapply(d,mean)
M
S <- cov(d)
S

# we build one of the possible contrast matrices to answer
# the question
C <- matrix(c(-1, 1, 0, 0,
              -1, 0, 1, 0,
              -1, 0, 0, 1), 3, 4, byrow=T)
C
# here we are looking at the effects on the d
# between the 1st and the 2nd, 3rd, 4th

# Test: H0: C%*%mu == 0 vs H1: C%*%mu != 0
alpha   <- .05
delta.0 <- c(0,0,0)

Md <- C %*% M 
Sd <- C %*% S %*% t(C)
Sdinv <- solve(Sd)

T2 <- n * t( Md - delta.0 ) %*% Sdinv %*% ( Md - delta.0 )

cfr.fisher <- ((q-1)*(n-1)/(n-(q-1)))*qf(1-alpha,(q-1),n-(q-1)) 

T2 < cfr.fisher
T2
cfr.fisher

# T2 is much higher than cfr.fisher => the p-value will be very small
P <- 1-pf(T2*(n-(q-1))/((q-1)*(n-1)),(q-1),n-(q-1))
P

### question (b)

# It is implicitly asking for confidence intervals on the components
# (for the mean of the increments after 8 hours, 16 hours and 24 hours)

# Simultaneous T2 intervals
IC.T2 <- cbind( Md - sqrt(cfr.fisher*diag(Sd)/n) , Md, Md + sqrt(cfr.fisher*diag(Sd)/n) )
IC.T2

# Bonferroni intervals 
k     <- q - 1   # number of increments (i.e., dim(C)[1])
cfr.t <- qt(1-alpha/(2*k),n-1)

IC.BF <- cbind( Md - cfr.t*sqrt(diag(Sd)/n) , Md, Md + cfr.t*sqrt(diag(Sd)/n) )
IC.BF


x11()
matplot(t(matrix(1:3,3,3)),t(IC.BF), type='b',pch='',xlim=c(0,4),xlab='',ylab='', main='Confidence intervals')
segments(matrix(1:3,3,1),IC.BF[,1],matrix(1:3,3,1),IC.BF[,3], col='orange', lwd=2)
points(1:3, IC.BF[,2], col='orange', pch=16)
points(1:3+.05, delta.0, col='black', pch=16)
segments(matrix(1:3+.1,3,1),IC.T2[,1],matrix(1:3+.1,3,1),IC.T2[,3], col='blue', lwd=2)
points(1:3+.1,IC.T2[,2], col='blue', pch=16)
legend('topright', c('Bonf. IC', 'Sim-T2 IC'), col=c('orange', 'blue'), lty=1, lwd=2)


### what happens if we change the constrast matrix?
Cbis <- matrix(c(-1, 1, 0, 0,
                 0, -1, 1, 0,
                 0, 0, -1, 1), 3, 4, byrow=T)
Cbis
# in this way we are looking at the mean increment of the d every 8 hours

Mdbis <- Cbis %*% M 
Sdbis <- Cbis %*% S %*% t(Cbis)
Sdinvbis <- solve(Sdbis)

T2bis <- n * t( Mdbis ) %*% Sdinvbis %*% Mdbis

T2bis < cfr.fisher

# compare the T2 test statistics associated with C and Cbis
T2bis
T2


# What is changed?
# The confidence intervals on the contrasts
# (because we are looking at different contrasts!)

IC.BFbis <- cbind( Mdbis - cfr.t*sqrt(diag(Sdbis)/n) , Mdbis, Mdbis + cfr.t*sqrt(diag(Sdbis)/n) )
IC.T2bis <- cbind( Mdbis - sqrt(cfr.fisher*diag(Sdbis)/n) , Mdbis, Mdbis + sqrt(cfr.fisher*diag(Sdbis)/n) )

IC.BFbis
IC.BF

IC.T2bis
IC.T2

### what if we want to verify the following hypothesis:
### "the drug decreases the d of two units with respect to
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

Md <- C %*% M 
Sd <- C %*% S %*% t(C)
Sdinv <- solve(Sd)

T2 <- n * t( Md - delta.0 ) %*% Sdinv %*% ( Md - delta.0 )

cfr.fisher <- ((q-1)*(n-1)/(n-(q-1)))*qf(1-alpha,(q-1),n-(q-1))
T2 < cfr.fisher
T2
cfr.fisher

# p-value
P <- 1-pf(T2*(n-(q-1))/((q-1)*(n-1)),(q-1),n-(q-1))
P


