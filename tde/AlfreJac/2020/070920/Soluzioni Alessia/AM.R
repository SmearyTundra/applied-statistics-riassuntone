########
# Problema 1
#######
data <- read.table('weather.txt',header = TRUE)
dim(data)
head(data)
cov(data)

data.std <- scale(data) #new dataset with mean zero and var=1
data.sd <- data.frame(data.std)

head(data.sd)

x11()
par(mar=rep(2,4))
boxplot(data.sd, las=2, col='gold')

pc.data <- princomp(data.sd, scores=T)
pc.data
summary(pc.data)


#Loadings
load.data <- pc.data$loadings
load.data

x11()
par(mar = c(2,2,2,1), mfrow=c(3,1))
for(i in 1:3) barplot(load.data[,i], ylim = c(-1, 1), main=paste('Loadings PC ',i,sep=''))

scores <- pc.data$scores #projection of the original dataset on the new reference system of PC
scores

#Representation of the scorses
x11()
par(mfrow=c(1,1))
plot(scores[,1:2])
abline(h=0, v=0, lty=2, col='grey') #each point correspond to a statistica unit and this is the projection
text(scores[,1:2])
head(data)

x11()
layout(matrix(c(2,3,1,3),2,byrow=T))
plot(pc.data, las=2, main='Principal Components', ylim=c(0,7)) #usefull to see the difference between the first two bar
abline(h=1, col='blue')
barplot(sapply(data.sd,sd)^2, las=2, main='Original Variables', ylim=c(0,7), ylab='Variances') #all equal to 1
plot(cumsum(pc.data$sde^2)/sum(pc.data$sde^2), type='b', axes=F, xlab='Number of components', ylab='Contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(data.sd),labels=1:ncol(data.sd),las=2)

summary(pc.data)


x_new <- data.frame(MeanTemp =30, MinTemp= 23, MaxTemp=36, DewPoint=22, Humidity= 65, Visibility=19,MeanWind= 5,MaxWind= 15)
x_new <- x_new - colMeans(data)/sapply(data,sd)
projection_1 <- as.numeric(x_new)%*%load.data[,1]
projection_2 <- as.numeric(x_new)%*%load.data[,2]
projection_3 <- as.numeric(x_new)%*%load.data[,3]
#######
# PROBLEMA 2
#####
data_c <- read.table('candle.txt')
head(data_c)
data_s <- read.table('sunshine.txt')
head(data_s)

load("~/Desktop/IV anno/Applied statistics/mcshapiro.test.RData")
mcshapiro.test(data_c)
mcshapiro.test(data_s)
cov(data_s)
cov(data_c)
t1 <- data_s
t2 <- data_c
n1 <- dim(t1)[1] 
n2 <- dim(t2)[1] 
p  <- dim(t1)[2] 


# we compute the sample mean, covariance matrices and the matrix Spooled
t1.mean <- sapply(t1,mean)
t2.mean <- sapply(t2,mean)
t1.cov  <-  cov(t1)
t2.cov  <-  cov(t2)
Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2) # S pooled
# we compare the matrices
list(S1=t1.cov, S2=t2.cov, Spooled=Sp)

Test
# Test H0: mu1 == mu2  vs  H1: mu1 != mu2
# i.e.,
# Test H0: mu1-mu2 == c(0,0)  vs  H1: mu1-mu2 != c(0,0)

alpha   <- .05
delta.0 <- c(0,0)
Spinv   <- solve(Sp)

T2 <- n1*n2/(n1+n2) * (t1.mean-t2.mean-delta.0) %*% Spinv %*% (t1.mean-t2.mean-delta.0)
#Malanobius distance induced by the S polled

cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
T2 < cfr.fisher

# p-value
P <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
P  

cfr.t <- qt(1-alpha/4,n1+n2-1-p)
IC.T2.X1 <- c(t1.mean[1]-t2.mean[1]-sqrt(cfr.t*Sp[1,1]*(1/n1+1/n2)), t1.mean[1]-t2.mean[1]+sqrt(cfr.t*Sp[1,1]*(1/n1+1/n2)) )
IC.T2.X2 <- c(t1.mean[2]-t2.mean[2]-sqrt(cfr.t*Sp[2,2]*(1/n1+1/n2)), t1.mean[2]-t2.mean[2]+sqrt(cfr.t*Sp[2,2]*(1/n1+1/n2)) )
IC.T2 <- rbind(IC.T2.X1, IC.T2.X2)
dimnames(IC.T2)[[2]] <- c('inf','sup')                        
IC.T2

a <- c(1,-1) 
cfr.t <- qt(1-alpha,n1+n2-1-p)
IC.T2.X1 <- c(t(a)%*%(t1.mean[1]-t2.mean[1]), t(a)%*%(t1.mean[1]-t2.mean[1])+sqrt(cfr.t*t(a)%*%Sp%*%a*(1/n1+1/n2)) )
dimnames(IC.T2.X1)[[2]] <- c('inf','sup')                        
IC.T2.X1

#####
#PROBLEMA 3
######
data <- read.table('leaven.txt')
head(data)
dim(data)
dummy <- rep(0,dim(data)[1])
attach(data)
yeast <- factor(yeast)
dummy[which(yeast=='sd')] <- 1

fm <- lm(volume ~ time +time: dummy +I(time^2) + dummy: I(time^2))
summary(fm)
coef <-fm$coefficients
coef[2]+coef[4]
coef[3]+coef[5]
x11()
par(mfrow=c(2,2))
plot(fm)
# If we have a path we need another regressor while if the sigma is differnt in different part we need to transform data. We can do the same with the studenties residuals ( we have also a red line related to the local mean trends ). Then we have QQ plot and the last plot is the lavarage related to the cook distance. We need to note the one that are out of the contour plot 
shapiro.test(residuals(fm))
library(MASS)
library(car)
library(rgl)

C <- c(0,0,1,0,0)
linearHypothesis(fm, C, 0)
C <- c(0,0,0,0,1)
linearHypothesis(fm, C, 0)

summary((fm))

fm <- lm(volume ~ time: dummy +I(time^2) + dummy: I(time^2))
summary(fm)
fm$coefficients[4]+ fm$coefficients[2]

Z0.new <- data.frame(time=2,dummy=0)
# new observation need to be in the shape of the dataframe with name coeherent with he one of the regression

# Conf. int. for the mean
Conf <- predict(fm, Z0.new, interval='confidence', level=1-0.01)  
Conf

#####
#PROBLEMA 4
####
library(fda)
noisycurve <- read.table("tide.txt",header=T)
head(noisycurve)
dim(noisycurve)
Xobs0 <- noisycurve$level
abscissa <- noisycurve$time
NT <- length(abscissa)

x11()
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
m <- 4
nbasis <- 6:30
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(min(abscissa),max(abscissa)), nbasis[i], m)
  gcv[i] <- smooth.basis(abscissa, Xobs0, basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis <- nbasis[which.min(gcv)]

help(create.bspline.basis) # we can have analogus function to get different basis
basis <- create.bspline.basis(rangeval=c(min(abscissa),max(abscissa)), nbasis=nbasis, norder=m)
# If breaks are not provided, equally spaced knots are created
names(basis)

plot(basis)

basismat <- eval.basis(abscissa, basis)
dim(basismat)
head(basismat)

help(lsfit)
lsfit(basismat, Xobs0, intercept=FALSE)$coef
Xsp0 <- basismat %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef
# final operation to do the fit 

x11()
par(mfrow=c(1,1))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)

S <- basismat%*%solve(t(basismat)%*%basismat)%*%t(basismat) #projection operator 
sigmahat <- sqrt(sum((Xsp0-Xobs0)^2)/(NT-nbasis)) #estimate of sigma
lb <- Xsp0-qnorm(0.975)*sigmahat*sqrt(diag(S))
ub <- Xsp0+qnorm(0.975)*sigmahat*sqrt(diag(S))

x11()
plot(abscissa,Xsp0,type="l",col="blue",lwd=2,ylab="")
points(abscissa,lb,type="l",col="blue",lty="dashed")
points(abscissa,ub,type="l",col="blue",lty="dashed")

# We want to keep under control also the derivatice of the data compute the central finite differences
rappincX1 <- (Xobs0[3:NT]-Xobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])
basismat1<- eval.basis(abscissa, basis, Lfdobj=1)
Xsp1 <- basismat1 %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef


plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
