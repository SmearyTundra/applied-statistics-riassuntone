#### FUNCTIONAL DATA ANALYSIS
library(fda)


########################### SMOOTHING ################################
#### ----------------------------------------
# Upload noisy data and visualise it
noisycurve <- read.table("noisycurvebis.txt",header=T)
head(noisycurve)
dim(noisycurve)

Xobs0 <- noisycurve$X0
abscissa <- noisycurve$Abscissa
NT <- length(abscissa) # number of locations of observations

x11()
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
plot(abscissa,Xobs0,xlab="t",ylab="observed data", type = "l")


# compute the central finite differences
rappincX1 <- (Xobs0[3:NT]-Xobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])
rappincX2 <- ((Xobs0[3:NT]-Xobs0[2:(NT-1)])/(abscissa[3:NT]-abscissa[2:(NT-1)])-(Xobs0[2:(NT-1)]-Xobs0[1:(NT-2)])/(abscissa[2:(NT-1)]-abscissa[1:(NT-2)]))*2/(abscissa[3:(NT)]-abscissa[1:(NT-2)])

par(mfrow=c(1,3))
plot(abscissa,Xobs0,xlab="t",ylab="observed data",type="l")
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")





#### REGRESSION SPLINES 
#### -------------------------------
# Set parameters
m <- 5           # spline order 
degree <- m-1    # spline degree 
nbasis <- 9

# Create the basis
basis <- create.bspline.basis(rangeval=range(abscissa), nbasis=nbasis, norder=m)
# If breaks are not provided, equally spaced knots are created
names(basis)
plot(basis)
# Evaluate the basis on the grid of abscissa
basismat <- eval.basis(abscissa, basis)
dim(basismat) # number of data x number of basis


# Fit via LS
est_coef = lsfit(basismat, Xobs0, intercept=FALSE)$coef

basismat <- eval.basis(abscissa, basis)
Xsp0 <- basismat %*% est_coef
# to obtain the first derivative (argument Lfdobj=1)
basismat1<- eval.basis(abscissa, basis, Lfdobj=1)
Xsp1 <- basismat1 %*% est_coef
# to obtain the second derivative (argument Lfdobj=2)
basismat2<- eval.basis(abscissa, basis, Lfdobj=2)
Xsp2 <- basismat2 %*% est_coef

par(mfrow=c(1,1))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
abline(v=basis$params)


x11(width = 14)
par(mfrow=c(1,3))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,truecurve$X0vera ,type="l",col="orange",lwd=3)
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
legend("topleft", legend = c("noisy data","true curve","estimated curve"), col = c("black", "orange","blue"), lwd = c(1,3,2))
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(truecurve$Abscissa,truecurve$X1vera,type='l',col="orange",lwd=3)
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(truecurve$Abscissa,truecurve$X2vera,type='l',col="orange",lwd=3)
points(abscissa,Xsp2 ,type="l",col="blue",lwd=2)


############################
# alternative code 
Xsp <- smooth.basis(argvals=abscissa, y=Xobs0, fdParobj=basis)
Xsp0bis <- eval.fd(abscissa, Xsp$fd) #  the curve smoothing the data
Xsp1bis <- eval.fd(abscissa, Xsp$fd, Lfd=1) # first derivative
Xsp2bis <- eval.fd(abscissa, Xsp$fd, Lfd=2) # second derivative
df <- Xsp$df   #  the degrees of freedom in the smoothing curve  
df             #  for regression splines the df are the number of basis
############################


#### Approximate pointwise confidence intervals 
# DA ESEGUIRE CON ALTERNATIVE CODE
# As in linear models, we can estimate the variance of x(t) as
# sigma^2*diag[phi*(phi'phi)^{-1}(phi)']
S <- basismat%*%solve(t(basismat)%*%basismat)%*%t(basismat) #projection operator 
sum(diag(S))
sigmahat <- sqrt(sum((Xsp0-Xobs0)^2)/(NT-df)) #estimate of sigma
lb <- Xsp0-qnorm(0.975)*sigmahat*sqrt(diag(S))
ub <- Xsp0+qnorm(0.975)*sigmahat*sqrt(diag(S))

x11()
plot(abscissa,Xsp0,type="l",col="blue",lwd=2,ylab="")
points(abscissa,lb,type="l",col="blue",lty="dashed")
points(abscissa,ub,type="l",col="blue",lty="dashed")
points(abscissa,truecurve$X0vera,type="l")



### generalized cross-validation for nbasis
nbasis <- 6:30
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(range(abscissa), nbasis[i], m)
  gcv[i] <- smooth.basis(abscissa, Xobs0, basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]
abline(v = nbasis[which.min(gcv)], col = 2)










#### FOURIER BASIS REGRESSION SPLINES
####---------------------------------------
# generalized cross-validation
nbasis <- 6:50
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.fourier.basis(range(abscissa), nbasis[i])
  gcv[i] <- smooth.basis(abscissa, Xobs0, basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]
abline(v=nbasis[which.min(gcv)],col='red')

basis <- create.fourier.basis(rangeval=range(abscissa), nbasis=nbasis[which.min(gcv)])
plot(basis)

Xsp <- smooth.basis(argvals=abscissa, y=Xobs0, fdParobj=basis)
Xsp0bis <- eval.fd(abscissa, Xsp$fd) #  the curve smoothing the data

plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0bis ,type="l",col="red",lwd=2)


# compute the central finite differences and first derivative
rappincX1 <- (Xobs0[3:NT]-Xobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])
Xsp1bis <- eval.fd(abscissa, Xsp$fd, Lfd=1) # first derivative

plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1bis,type='l',col="orange",lwd=3)



### Overfitting
nbasis <- 50
basis <- create.fourier.basis(rangeval=range(abscissa), nbasis=nbasis)

Xsp <- smooth.basis(argvals=abscissa, y=Xobs0, fdParobj=basis)
Xsp0bis <- eval.fd(abscissa, Xsp$fd) #  the curve smoothing the data

plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0bis ,type="l",col="red",lwd=2)




### Oversmoothing
nbasis <- 5
basis <- create.fourier.basis(rangeval=range(abscissa), nbasis=nbasis)

Xsp <- smooth.basis(argvals=abscissa, y=Xobs0, fdParobj=basis)
Xsp0bis <- eval.fd(abscissa, Xsp$fd) #  the curve smoothing the data

plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0bis ,type="l",col="red",lwd=2)






















########################### SMOOTHING SPLINES WIHT FDPAR and LAMBDA ################################
#### ----------------------------------------
# riga 270 FDA_SMOOTHING
# breaks <- abscissa[((0:50)*2)+1]
breaks <- abscissa

basis <- create.bspline.basis(breaks, norder=m)
functionalPar <- fdPar(fdobj=basis, Lfdobj=3, lambda=1e-8)  
# functional parameter, having arguments: 
# basis, order of the derivative to be penalized, smoothing parameter.

Xss <- smooth.basis(abscissa, Xobs0, functionalPar)

Xss0 <- eval.fd(abscissa, Xss$fd, Lfd=0)
Xss1 <- eval.fd(abscissa, Xss$fd, Lfd=1)
Xss2 <- eval.fd(abscissa, Xss$fd, Lfd=2)

df <- Xss$df   #  the degrees of freedom in the smoothing curve
df
gcv <- Xss$gcv  #  the value of the gcv statistic
gcv

x11()
par(mfrow=c(1,3))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xss0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xss1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xss2 ,type="l",col="blue",lwd=2)

##### Using cross-validation

lambda <- 10^seq(-12,-5,by = 0.5) # grid
gcv <- numeric(length(lambda))
for (i in 1:length(lambda)){
  functionalPar <- fdPar(fdobj=basis, Lfdobj=3, lambda=lambda[i])  
  gcv[i] <- smooth.basis(abscissa, Xobs0, functionalPar)$gcv
}
par(mfrow=c(1,1))
plot(log10(lambda),gcv)
lambda[which.min(gcv)]


# best lambda
functionalParbest <- fdPar(fdobj=basis, Lfdobj=3, lambda=lambda[which.min(gcv)])  

Xssbest <- smooth.basis(abscissa, Xobs0, functionalParbest)

Xss0best <- eval.fd(abscissa, Xssbest$fd, Lfd=0)
Xss1best <- eval.fd(abscissa, Xssbest$fd, Lfd=1)
Xss2best <- eval.fd(abscissa, Xssbest$fd, Lfd=2)

dfbest <- Xssbest$df   #  the degrees of freedom in the smoothing curve
dfbest
gcvbest <- Xssbest$gcv  #  the value of the gcv statistic
gcvbest


par(mfrow=c(1,3))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xss0ter ,type="l",col="red",lwd=1)
points(abscissa,Xss0bis ,type="l",col="green",lwd=1)
points(abscissa,Xss0best ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xss1ter ,type="l",col="red",lwd=1)
points(abscissa,Xss1bis ,type="l",col="green",lwd=1)
points(abscissa,Xss1best ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xss2ter ,type="l",col="red",lwd=1)
points(abscissa,Xss2bis ,type="l",col="green",lwd=1)
points(abscissa,Xss2best ,type="l",col="blue",lwd=2)








########################### LOCAL POLY REGRESSION ################################
#### ----------------------------------------
#### riga 389 FDA_SMOOTHING

library(KernSmooth)

# Xobs0 era una singola curva
m <- 5           # order of the polynomial
degree <- m-1    # degree of the polynomial


bw <- 0.05 # bandwidth
# a too large bandwidth leads to oversmoothing
# a too small bandwidth leads to undersmoothing
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Recommendation: when choosing the bandwidth, look at the
# derivatives vs central finite differences



Xsm0 <- locpoly(abscissa, Xobs0, degree=degree,
                bandwidth=bw, gridsize=length(abscissa), 
                range.x=range(abscissa))
Xsm0 <- Xsm0$y

par(mfrow=c(1,1))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsm0 ,type="l",col="blue")

# Estimate of the derivatives

Xsm1 <- locpoly(abscissa,Xobs0,drv=1,degree=degree,bandwidth=bw,
                gridsize=length(abscissa), range.x=range(abscissa))
Xsm1 <- Xsm1$y

Xsm2 <- locpoly(abscissa,Xobs0,drv=2,degree=degree,bandwidth=bw,
                gridsize=length(abscissa), range.x=range(abscissa))
Xsm2 <- Xsm2$y


par(mfrow=c(1,3))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsm0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsm1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsm2 ,type="l",col="blue",lwd=2)









########################### COSTRAINED FUNCTIONS ################################
#### ----------------------------------------
#### riga 493














