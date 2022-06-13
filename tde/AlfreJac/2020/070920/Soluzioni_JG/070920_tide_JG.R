tide <- read.table("tide.txt", header = T)
library(fda)
m <- 4           # spline order 
degree <- m-1    # spline degree 

abscissa <- tide$time
Xobs0 <- tide$level


nbasis <- 4:45
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(0,23.5), nbasis[i], m)
  gcv[i] <- smooth.basis(abscissa, Xobs0, basis)$gcv
}
x11()
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]
#nbasis = 12

basis <- create.bspline.basis(c(0,23.5), nbasis[which.min(gcv)], m)
basismat <- eval.basis(abscissa, basis)

lsfit(basismat, Xobs0, intercept=FALSE)$coef #Intercept = F since basis already generate it

Xsp0 <- basismat %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef

x11()
plot(basis)
x11()
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)

#we consider smoothing as well
lambda <- c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12)
gcv <- numeric(length(lambda))

for (i in 1:length(lambda)){
  functionalPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda[i])  
  gcv[i] <- smooth.basis(abscissa, Xobs0, functionalPar)$gcv
}
x11()
par(mfrow=c(1,1))
plot(log10(lambda),gcv)
lambda[which.min(gcv)]

Xss <- smooth.basis(abscissa, Xobs0, fdPar(fdobj=basis, Lfdobj=2, lambda=lambda[which.min(gcv)]))

Xss0 <- eval.fd(abscissa, Xss$fd, Lfd=0)
points(abscissa,Xss0 ,type="l",col="red",lwd=2)
#with a penalization factor we do not change the smoothing that much


S <- basismat%*%solve(t(basismat)%*%basismat)%*%t(basismat) #projection operator 
sigmahat <- sqrt(sum((Xsp0-Xobs0)^2)/(NT-nbasis)) #estimate of sigma
lb <- Xsp0-qnorm(0.975)*sigmahat*sqrt(diag(S))
ub <- Xsp0+qnorm(0.975)*sigmahat*sqrt(diag(S))

x11()
plot(abscissa,Xsp0,type="l",col="blue",lwd=2,ylab="")
points(abscissa,lb,type="l",col="blue",lty="dashed")
points(abscissa,ub,type="l",col="blue",lty="dashed")



rappincX1 <- (Xobs0[3:NT]-Xobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])  # [y(t+1) - y(t-1)]/2*deltat =~ y'(t)

basismat1<- eval.basis(abscissa, basis, Lfdobj=1)
Xsp1 <- basismat1 %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef
x11()
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
#of course the finite difference approximation is very rough, but still the fit is acceptable
#we see however that the derivative of our fit is not as smooth as one should expect it to be
#an appropriate basis for this phenomenon coul be a fourier basis, given the periodic nature of the tide
#such basis moreover preserves the oscillatory behaviour in the derivative as well