#######
# MULTIPLE POPULATIONS
######

library(car)

###-----------------------------------------------------------------------------------------
##### INDEPENDENT - Test on mean, Bonferroni IC, pvalue, ellipse, test on reductions #####
###-----------------------------------------------------------------------------------------
rm(list=ls())
load('mcshapiro.test.RData')
tilda = '~'
graphics.off()

candle <- read.table('candle.txt', header = TRUE)
sunshine <- read.table('sunshine.txt', header = TRUE)
head(candle)
head(sunshine)

mcshapiro.test(candle)$pvalue
mcshapiro.test(sunshine)$pvalue

n1 <- dim(candle)[1]
n2 <- dim(sunshine)[1]
p <- dim(candle)[2]
alpha <- 0.05

# Test H0: mu1-mu2 == c(0,0)  vs  H1: mu1-mu2 != c(0,0)

m1 <- colMeans(candle)
m2 <- colMeans(sunshine)

S1 <- cov(candle)
S2 <- cov(sunshine)
Sp <- (S1 * (n1 - 1) + S2 * (n2 - 1)) / (n1 + n2 -2)
invSP <- solve(Sp)

T2 <- n1*n2/(n1 + n2)*(m1 - m2)%*%invSP%*%(m1 - m2)
cfr.fisher <- p*(n1 + n2 - 2)/(n1 + n2 - 1 - p)*qf(1-alpha, p, n1 + n2 - 1 - p)

T2 < cfr.fisher

### Ellipsoidal region
plot(candle-sunshine, asp = 1)
ellipse(m1 - m2, shape=Sp/(n1 + n2), sqrt(cfr.fisher), col = 'blue', lty = 2, center.pch = 16)
points(0, 0, pch = 16, col = 'red', cex=1.5)


# Compute the p-value
P <- 1-pf(T2*(n1 + n2 - 1 - p)/(p*(n1 + n2 - 2)), p, n1 + n2 - 1 - p)
P

#Bonferroni intervals
k <- 2
cfr.t <- qt(1-alpha/(2*k), n1 + n2 - 2)
IC1 <- c( m1[1] - m2[1] - sqrt(Sp[1,1] * (n1 + n2)/(n1*n2)) * cfr.t, m1[1] - m2[1],  m1[1] - m2[1] + sqrt(Sp[1,1] * (n1 + n2)/(n1*n2)) * cfr.t )
IC2 <- c( m1[2] - m2[2] - sqrt(Sp[2,2] * (n1 + n2)/(n1*n2)) * cfr.t, m1[2] - m2[2],  m1[2] - m2[2] + sqrt(Sp[2,2] * (n1 + n2)/(n1*n2)) * cfr.t )

#reduction
candle.red <- candle[,1] - candle[,2]
sunshine.red <- sunshine[,1] - sunshine[,2]
d1 <- mean(candle.red)
d2 <- mean(sunshine.red)

t.test(candle.red, sunshine.red, mu = 0, alternative = 'greater', conf.level = 0.95)



###-------------------------------------------------------------------------
#### PAIRED - Test on mean, simultaneous T2 intervals, Bonferroni IC, pvalue 
###------------------------------------------------------------------
rm(list=ls())
load('mcshapiro.test.RData')
tilda = '~'
graphics.off()

effluent <- read.table('effluent.dat', header=T)
head(effluent)

colnames(effluent) <- c('BOD_Lab1','SS_Lab1','BOD_Lab2','SS_Lab2')

pairs(effluent,pch=19, main='Dataset effluent')

dev.off()

# we compute the sample of differences
D <- data.frame(DBOD=effluent[,1]-effluent[,3], DSS=effluent[,2]-effluent[,4]) 
D

plot(D, asp=1, pch=19, main='Dataset of Differences')
abline(h=0, v=0, col='grey35')
points(0,0, pch=19, col='grey35')

dev.off()

### T2 Hotelling Test 
# H0: delta == delta.0 vs H1: delta != delta.0 with delta.0=c(0,0)

# Test the Gaussian assumption (on D!)
mcshapiro.test(D)

n <- dim(D)[1]  
p <- dim(D)[2]  

D.mean   <- sapply(D,mean)
D.cov    <- cov(D)
D.invcov <- solve(D.cov)

alpha   <- .05
delta.0 <- c(0,0)

D.T2 <- n * (D.mean-delta.0) %*% D.invcov %*% (D.mean-delta.0)
D.T2

cfr.fisher <- (n-1)*p/(n-p)*qf(1-alpha,p,n-p)
cfr.fisher

D.T2 < cfr.fisher # if FALSE we reject H0 at level 5%

# we compute the p-value
P <- 1-pf(D.T2*(n-p)/(p*(n-1)), p, n-p)
P

# Ellipsoidal confidence region with confidence level 95%
plot(D, asp=1, pch=1, main='Dataset of the Differences',ylim=c(-15,60))
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt(cfr.fisher), lwd=2)

### Simultanous T2 intervals
IC.T2.DBOD <- c( D.mean[1]-sqrt(cfr.fisher*D.cov[1,1]/n) , D.mean[1], D.mean[1]+sqrt(cfr.fisher*D.cov[1,1]/n) )
IC.T2.DSS  <- c( D.mean[2]-sqrt(cfr.fisher*D.cov[2,2]/n) , D.mean[2], D.mean[2]+sqrt(cfr.fisher*D.cov[2,2]/n) )

T2 <- rbind(IC.T2.DBOD, IC.T2.DSS)
dimnames(T2)[[2]] <- c('inf','center','sup')
T2

### Bonferroni intervals
k <- p  # 2
cfr.t <- qt(1-alpha/(2*k),n-1)

IC.BF.DBOD <- c( D.mean[1]-cfr.t*sqrt(D.cov[1,1]/n) , D.mean[1], D.mean[1]+cfr.t*sqrt(D.cov[1,1]/n) )
IC.BF.DSS  <- c( D.mean[2]-cfr.t*sqrt(D.cov[2,2]/n) , D.mean[2], D.mean[2]+cfr.t*sqrt(D.cov[2,2]/n) )

Bf <- rbind(IC.BF.DBOD, IC.BF.DSS)
dimnames(Bf)[[2]] <- c('inf','center','sup')
Bf

-----
###---------------------------------------------------------------------------
#### REPEATED - Test on mean, simultaneous T2 intervals, Bonferroni IC, pvalue
###----------------------------------------------------------------------------

pressure <- read.table ('pressure.txt', col.names=c('h.0','h.8','h.16','h.24'), header=T)
head(pressure)
dim(pressure)

mcshapiro.test(pressure)

matplot(t(pressure), type='l')

# (a) Perform a test at level 5% to prove that the drug has influence on 
#     the blood pressure during the 24 hours

n <- dim(pressure)[1]
p <- dim(pressure)[2]

M <- sapply(pressure,mean)
S <- cov(pressure)

# we build one of the possible contrast matrices to answer the question
C <- matrix(c(-1, 1, 0, 0,
              -1, 0, 1, 0,
              -1, 0, 0, 1), 3, 4, byrow=T)

# Test: H0: C%*%mu == 0 vs H1: C%*%mu != 0
alpha   <- .05
delta.0 <- c(0,0,0)

Md <- C %*% M 
Sd <- C %*% S %*% t(C)
Sdinv <- solve(Sd)

T2 <- n * t( Md - delta.0 ) %*% Sdinv %*% ( Md - delta.0 )

cfr.fisher <- ((p-1)*(n-1)/(n-(p-1)))*qf(1-alpha,(p-1),n-(p-1)) 

T2 < cfr.fisher

P <- 1-pf(T2*(n-(p-1))/((p-1)*(n-1)),(p-1),n-(p-1))
P

# (b) Highlight the effect of the drug on the blood pressure
# It is implicitly asking for confidence intervals on the components
# (for the mean of the increments after 8 hours, 16 hours and 24 hours)

# Simultaneous T2 intervals
IC.T2 <- cbind( Md - sqrt(cfr.fisher*diag(Sd)/n) , Md, Md + sqrt(cfr.fisher*diag(Sd)/n) )
IC.T2

# Bonferroni intervals 
k     <- p - 1   # number of increments (i.e., dim(C)[1])
cfr.t <- qt(1-alpha/(2*k),n-1)

IC.BF <- cbind( Md - cfr.t*sqrt(diag(Sd)/n) , Md, Md + cfr.t*sqrt(diag(Sd)/n) )
IC.BF

#Bonferroni for variance 
ICBVar<-data.frame('Inf'=(n-1)*diag(S)/qchisq(1-alpha/(2*k),n-1),
                   'M'  =diag(S),
                   'Sup'=(n-1)*diag(S)/qchisq(alpha/(2*k),n-1))

#########
### what if we want to verify the following hypothesis:
### "the drug decreases the pressure of two units with respect to
### the baseline at both 8 and 16 hours, and its effect vanishes in 24 hours
### from the drug administration"

C <- matrix(c(-1, 1, 0, 0,
              -1, 0, 1, 0,
              -1, 0, 0, 1), 3, 4, byrow=T)
delta.0 <- c(-2,-2,0)




