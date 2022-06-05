#########
# FUNCTIONAL
#########

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

library(fdakma)
fdakma_example_noalign_0der <- kma(
  x=1:80, y0=t(spectra), n.clust = 3, 
  warping.method = 'affine', 
  similarity.method = 'd0.pearson',   # similarity computed as the cosine
  # between the original curves 
  # (correlation)
  center.method = 'k-means'
  #,seeds = c(1,11,21) # you can give a little help to the algorithm...
)

kma.show.results(fdakma_example_noalign_0der)
fdakma_example_noalign_0der$labels

