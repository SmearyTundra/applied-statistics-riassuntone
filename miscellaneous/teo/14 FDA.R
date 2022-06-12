#########
# FUNCTIONAL
#########

#### Reggression splines
library(fda)


norder <- 5      # spline order (4th order polynomials)
degree <- m-1    # spline degree
nbasis <- 9      # how many basis we want

basis <- create.bspline.basis(rangeval=c(0, 1),
                              nbasis=nbasis,
                              norder=norder)

basismat <- eval.basis(abscissa, basis)
est_coef <- lsfit(basismat, Xobs0, intercept=FALSE)$coef
Xsp0 <- basismat %*% est_coef

# plot fitted line
plot(abscissa, Xobs0, xlab="t", ylab="observed data")
points(abscissa, Xsp0, type="l", col="blue", lwd=2)
abline(v=basis$params)
# 1st derivative
basismat1 <- eval.basis(abscissa, basis, Lfdobj=1)
Xsp1 <- basismat1 %*% est_coef
# 2nd derivative
basismat2 <- eval.basis(abscissa, basis, Lfdobj=2)
Xsp2 <- basismat2 %*% est_coef


# faster
Xsp <- smooth.basis(argvals=abscissa, y=Xobs0, fdParobj=basis)
Xsp0bis <- eval.fd(abscissa, Xsp$fd) #  the curve smoothing the data
Xsp1bis <- eval.fd(abscissa, Xsp$fd, Lfd=1) # first derivative
Xsp2bis <- eval.fd(abscissa, Xsp$fd, Lfd=2) # second derivative



# basis choice with cross validation
nbasis <- 6:30
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
    basis <- create.bspline.basis(c(0, 1), nbasis[i], m)
    gcv[i] <- smooth.basis(abscissa, Xobs0, basis)$gcv
}
par(mfrow=c(1, 1))
plot(nbasis, gcv)
nbasis[which.min(gcv)]
abline(v=nbasis[which.min(gcv)], col=2)


# specify breaks
breaks <- abscissa
basis <- create.bspline.basis(breaks, norder=m)
functionalPar <- fdPar(fdobj=basis, Lfdobj=3, lambda=1e-8) # try changing lambda
                                                            # 1e-5 oversmoothing
                                                            # 1e-12 overfit
Xss <- smooth.basis(abscissa, Xobs0, functionalPar)
Xss0 <- eval.fd(abscissa, Xss$fd, Lfd=0)
Xss1 <- eval.fd(abscissa, Xss$fd, Lfd=1)
Xss2 <- eval.fd(abscissa, Xss$fd, Lfd=2)

# and use CV to set lambda
lambda <- 10^seq(-12, -5, by=0.5)
gcv <- numeric(length(lambda))
for (i in 1:length(lambda)){
    functionalPar <- fdPar(fdobj=basis, Lfdobj=3, lambda=lambda[i])
    gcv[i] <- smooth.basis(abscissa, Xobs0, functionalPar)$gcv
}
par(mfrow=c(1, 1))
plot(log10(lambda), gcv)
lambda[which.min(gcv)]

# substitute it
functionalParbest <- fdPar(fdobj=basis, Lfdobj=3, lambda=lambda[which.min(gcv)])
Xssbest <- smooth.basis(abscissa, Xobs0, functionalParbest)
Xss0best <- eval.fd(abscissa, Xssbest$fd, Lfd=0)
Xss1best <- eval.fd(abscissa, Xssbest$fd, Lfd=1)
Xss2best <- eval.fd(abscissa, Xssbest$fd, Lfd=2)




#### Local Polynomial Regression
library(KernSmooth)
# in lab 14 copia tutto












##### FPCA - BSPLINE + HARMONICS + SCORES  ####
rm(list = ls())
graphics.off()
traffic <- read.table('traffic.txt')
traffic <- t(traffic)
matplot(traffic, type = 'l')

nbasis <- 15
time <- 1:24
basis <- create.bspline.basis(rangeval=c(1,24), nbasis=nbasis, norder=4)
data_W.fd.1 <- Data2fd(y = traffic,argvals = time,basisobj = basis)
plot.fd(data_W.fd.1)

data_W.fd.1$coefs

#FPCA
pca_W.1 <- pca.fd(data_W.fd.1,nharm=5,centerfns=TRUE)
pca_W.1$varprop
# scree plot
# pca.fd computes all the 365 eigenvalues, but only the first N-1=34 are non-null
plot(pca_W.1$values[1:15],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_W.1$values)[1:15]/sum(pca_W.1$values),xlab='j',ylab='CPV',ylim=c(0.8,1))

par(mfrow=c(1,3))
plot(pca_W.1$harmonics[1,],col=1,ylab='FPC1',ylim=c(-0.3,0.4))
abline(h=0,lty=2)
plot(pca_W.1$harmonics[2,],col=2,ylab='FPC2',ylim=c(-0.3,0.4))
plot(pca_W.1$harmonics[3,],col=1,ylab='FPC3',ylim=c(-0.3,0.4))

par(mfrow=c(1,3))
plot.pca.fd(pca_W.1, nx=100, pointplot=TRUE, harm=c(1,2,3), expand=0, cycle=FALSE)
dev.off()

plot(pca_W.1$scores[,1], pch = 16, xlab = 'DAY', ylab = 'PC1')
plot(pca_W.1$scores[,2], pch = 16, xlab = 'DAY', ylab = 'PC2')

plot(pca_W.1$scores[,1],pca_W.1$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)

##### CLUSTERING ####
rm(list=ls())
spectra <- read.table('spectra.txt')
spectra <- t(spectra)
x11()
matplot(spectra, type = 'l')

nbasis <- 5:30
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(1,80), nbasis[i], norder = 4)
  gcv[i] <- smooth.basis(1:80, spectra[,1], basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]
basis <- create.bspline.basis(c(1,80), 11, norder = 4)
#mode 1
basismat <- eval.basis(1:80, basis)
lsfit(basismat, spectra[,1], intercept=FALSE)$coef
#mode2
coefs <- smooth.basis(1:80, spectra[,1], basis)$fd$coefs[1:3]

wl <- 1:80
data_W.fd.1 <- Data2fd(y = spectra,argvals = wl,basisobj = basis)
plot.fd(data_W.fd.1)




# kmean alignment

library(fdakma)
fdakma_example <- kma(
  x=1:80, y0=t(spectra), n.clust = 3, 
  warping.method = 'NOalignment-affine', 
  similarity.method = 'd0.pearson-d1.pearson',   # similarity computed as the cosine
  # between the original curves 
  # (correlation)
  # similarity related to which derivative
  center.method = 'k-means'
  #,seeds = c(1,11,21) # you can give a little help to the algorithm...
)

kma.show.results(fdakma_example)
fdakma_example$labels

#fdakma_example$shift
#fdakma_example$dilation

# compare different warping and cluster
# plotta in automatico
kma.compare_example_3 <- kma.compare (
 x=x, y0=y0, y1=y1, n.clust=1:3,
 warping.method=c("NOalignment", "shift", "dilation", "affine"),
 similarity.method='d1.pearson',
 center.method='k-means',
 seeds=c(1, 21, 30),
 plot.graph=TRUE)