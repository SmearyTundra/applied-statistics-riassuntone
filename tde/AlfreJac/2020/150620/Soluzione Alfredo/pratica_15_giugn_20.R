#########33
## Problem 1.
##########
setwd("Documents/Polimi/APPSTATS")

pol <- read.table("pollution.txt")
head(pol)
n <- dim(pol)[1]
p <- dim(pol)[2]
n
p
mcshapiro.test(pol)
# assumption of normality verified. 
mu0 <- c(50,50)
X <- sapply(pol, mean)
S <- cov(pol)
T2 <- n * t(X-mu0) %*% solve(S) %*% (X-mu0)
alpha = 0.05
cfr.F <-  ( (n-1)*p / (n-p) )* qf(1-alpha, p, n-p)
T2 < cfr.F
P <- 1-pf(T2*(n-p)/(p*(n-1)), p, n-p)
P

# b)
library(car)
x11()
plot(pol, asp=1, ylim=c(-50,450),pch=1, main='Pollution')
ellipse(center=X, shape=S, radius=sqrt(cfr.F/n), lwd=2)
points(X[1],X[2], pch=16, col='grey35', cex=1.5)    
points(mu0[1],mu0[2], pch=16, col='pink1', cex=1.5)   
2*sqrt(cfr.F/n)*sqrt(eigen(S)$values) 
eigen(S)$vectors

# c
"Falls outside rejection region, that's why we reject."


# d)
IC.T2.X1 <- c(X[1]-sqrt(cfr.F*S[1,1]/n), X[1]+sqrt(cfr.F*S[1,1]/n), X[1]) 
IC.T2.X2 <- c(X[2]-sqrt(cfr.F*S[2,2]/n), X[2]+sqrt(cfr.F*S[2,2]/n), X[2]) 
IC.T2.X1
IC.T2.X2



### 
" Problem 2"
stone <- read.table("stoneflakes.txt")
head(stone)
## a)
st.e <-  dist(stone, method='euclidean')
st.ew <- hclust(st.e, method='ward.D2')
plot(st.ew, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(st.ew, k=3)

# b)
cluster.ec <- cutree(st.ew, k=3)
cluster.ec
# homoschedasticity
library(biotools)
boxM(stone, cluster.ec)
# normality of each group
mcshapiro.test(stone[which(cluster.ec==1),])
mcshapiro.test(stone[which(cluster.ec==2),])
mcshapiro.test(stone[which(cluster.ec==3),])
n1 <- dim(stone[which(cluster.ec==1),])[1]
n2 <- dim(stone[which(cluster.ec==2),])[1]
n3 <- dim(stone[which(cluster.ec==3),])[1]


fit <- manova(as.matrix(stone) ~ as.factor(cluster.ec))
summary.manova(fit,test="Hotelling-Lawley")
shapiro.test(fit$residuals)


# c)
alpha <- 0.1
g <- 3
k <- p*g*(g-1)/2
qT <- qt(1-alpha/(2*k), n-g)
m  <- sapply(stone,mean)         # estimates mu
m1 <- sapply(stone[which(cluster.ec==1),],mean)    # estimates mu.1=mu+tau.1
m2 <- sapply(stone[which(cluster.ec==2),],mean)    # estimates mu.2=mu+tau.2
m3 <- sapply(stone[which(cluster.ec==3),],mean)    # estimates mu.3=mu+tau.3  
W <- t(fit$residuals)%*%(fit$residuals)
inf12 <- m1-m2 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
sup12 <- m1-m2 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
inf13 <- m1-m3 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n3) )
sup13 <- m1-m3 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n3) )
inf23 <- m2-m3 - qT * sqrt( diag(W)/(n-g) * (1/n2+1/n3) )
sup23 <- m2-m3 + qT * sqrt( diag(W)/(n-g) * (1/n2+1/n3) )
CI <- list(C1_C2=cbind(inf12, sup12), C1_C3=cbind(inf13, sup13), C2_C3=cbind(inf23, sup23))
CI


### PRoblem 3
air <- read.table("airfoil.txt")
air$velocity <- as.factor(air$velocity)
head(air)
lm1 <- lm(sound ~velocity + velocity:frequency, data=air)
summary(lm1)
plot(lm1)
shapiro.test(lm1$residuals)
mcshapiro.test(lm1$residuals)

# statistical tests:
library(car)
C <- c(0,0,1,0)
linearHypothesis(lm1, C)

C <- matrix(c(0,1,0,0,
              0,0,1,0,
              0,0,0,1,
              1,0,0,0), byrow = T, nrow = 4, ncol = 4)
linearHypothesis(lm1, C)

C <- c(0,0,1,-1)
linearHypothesis(lm1, C)
"Not significant"

#c)
lm2 <- lm(lm(sound ~velocity + frequency, data=air))
summary(lm2)
lm2$coefficients

#d)
z0 <- data.frame(velocity="H", frequency=15000)
predict(lm2, z0, interval='confidence', level=1-0.05)  


###################
### Problem 4
################
rev <- read.table("revenues.txt")
head(rev)
attach(rev)
library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics
coordinates(rev) <- c('x','y')

## a)
# for the GLS we need a model for the variogram (i.e. for the residuals)
rvgm <- variogram(revenue ~ 1, rev)
plot(rvgm)
fit.r <- fit.variogram(rvgm, vgm(800, "Sph", range=2000))
plot(rvgm, fit.r)

# fit gls
gls.r <- gstat(formula = revenue ~ 1+population,
               data = rev,nmax=100, model=fit.r,
               set = list(gls=1))
summary(gls.r)
# Model assumptions: E[Zs] = ms (stationarity)
# It minimizes the mahalanobis distance induced by covariance of residuals
# covariance of residuals estimated by variogram, spherical model no nugget


## b)
# linear model to estimate population
lmaux <- lm(population~distance)
summary(lmaux)
shapiro.test(lmaux$residuals)
# assume linear model independent errors
# the errors are normally distributed so assumption is complied with
diff <- c(514711.6, 5033903.0) - c(514703.8,5035569.3)
z0 <- data.frame(distance=sqrt(t(diff) %*% diff))
predict(lmaux, z0, interval='prediction', level=1-0.05)  

zs0 <- data.frame(x=514703.8, y=5035569.3, population=6132.345)
coordinates(zs0) <- c('x','y')
predict(gls.r, zs0)

# c) kriging variance is 199.3191. It is not,as we also have uncertainty
# of the linear model to obtian the population and of the estimator of the
# residuals delta(si) covariance (variogram model)

