#################### FUNCTIONAL PCA
###----------------------------------- 

library(fda)
library(fields)



#### First: import and smoothing data
###------------------------------------
data_W <- CanadianWeather$dailyAv[,,1]
head(data_W)
dim(data_W)
matplot(data_W,type='l',main='Canadian temperature',xlab='Day',ylab='Temperature')

# First of all we smooth the data. We choose a Fourier basis
# (periodic). We need to set the dimension of the basis
time <- 1:365
# Choice 1: we set a high dimensional basis (interpolating)
# Pros: no loss of information
# Cons: possible overfitting 
basis.1 <- create.fourier.basis(rangeval=c(0,365),nbasis=365)
data_W.fd.1 <- Data2fd(y = data_W,argvals = time,basisobj = basis.1)
plot.fd(data_W.fd.1)

# Choice 2: reduced dimensionality (we set a low dimensional basis)
# Pros: the data are much smoother and the measurement error is filtered
# Cons: I could have lost important information
basis.2 <- create.fourier.basis(rangeval=c(0,365),nbasis=21)
data_W.fd.2 <- Data2fd(y = data_W,argvals = time,basisobj = basis.2)
plot.fd(data_W.fd.2)

# Choice 3: compromise between 1 and 2
basis.3 <- create.fourier.basis(rangeval=c(0,365),nbasis=109)
data_W.fd.3 <- Data2fd(y = data_W,argvals = time,basisobj = basis.3)
plot.fd(data_W.fd.3)





##### PERFORMING FPCA
###-------------------------------
# interpolated data (Choice 1)
plot.fd(data_W.fd.1,ylab='temperature')

pca_W.1 <- pca.fd(data_W.fd.1,nharm=5,centerfns=TRUE)

# scree plot
# pca.fd computes all the 365 eigenvalues, but only the first 
# N-1=34 are non-null
plot(pca_W.1$values[1:35],xlab='j',ylab='Eigenvalues') ### CHANGE 1:nsplines
plot(cumsum(pca_W.1$values)[1:35]/sum(pca_W.1$values),xlab='j',ylab='CPV',ylim=c(0.8,1))

# first two FPCs
x11()
par(mfrow = c(1,2))
plot(pca_W.1$harmonics[1,],col=1,ylab='FPC1')
abline(h=0,lty=2)
plot(pca_W.1$harmonics[2,],col=2,ylab='FPC2')

# plot of the FPCs as perturbation of the mean
media <- mean.fd(data_W.fd.1)

plot(media,lwd=2,ylim=c(-25,20),ylab='temperature',main='FPC1')
lines(media+pca_W.1$harmonics[1,]*sqrt(pca_W.1$values[1]), col=2)
lines(media-pca_W.1$harmonics[1,]*sqrt(pca_W.1$values[1]), col=3)
# variation in amplitude (more in winter than in summer)

plot(media,lwd=2,ylim=c(-20,20),ylab='temperature',main='FPC2')
lines(media+pca_W.1$harmonics[2,]*sqrt(pca_W.1$values[2]), col=2)
lines(media-pca_W.1$harmonics[2,]*sqrt(pca_W.1$values[2]), col=3)
# temperate climate or not

# Command of the library fda that automatically does these plots
par(mfrow=c(1,2))
plot(pca_W.1, nx=100, pointplot=TRUE, harm=c(1,2), expand=0, cycle=FALSE)


# scatter plot of the scores
par(mfrow=c(1,2))
plot(pca_W.1$scores[,1],pca_W.1$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
points(pca_W.1$scores[35,1],pca_W.1$scores[35,2],col=2, lwd=4) # higlightin a particular object


plot(pca_W.1$scores[,1],pca_W.1$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2",xlim=c(-400,250))
text(pca_W.1$scores[,1],pca_W.1$scores[,2],dimnames(data_W)[[2]], cex=1)

