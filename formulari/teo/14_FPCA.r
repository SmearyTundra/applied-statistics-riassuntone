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
# or t(data_W)

#### BSPLINE
# set the parameters
norder <- 5         # spline order (4th order polynomials)
degree <- norder-1  # spline degree
nbasis <- 9         # how many basis we want

time <- 1:24
abscissa <- time

basis <- create.bspline.basis(rangeval=range(abscissa), # 
                              nbasis=nbasis,
                              norder=norder)
data_W.bspline <- Data2fd(y = data_W,argvals = time,basisobj = basis)
plot.fd(data_W.bspline)


#### FOURIER

# First of all we smooth the data. We choose a Fourier basis
# (periodic). We need to set the dimension of the basis
time <- 1:365
abscissa <- time
# Choice 1: we set a high dimensional basis (interpolating)
# Pros: no loss of information
# Cons: possible overfitting 
basis.1 <- create.fourier.basis(rangeval=range(abscissa),nbasis=365)
data_W.fd.1 <- Data2fd(y = data_W,argvals = time,basisobj = basis.1)
plot.fd(data_W.fd.1)

# Choice 2: reduced dimensionality (we set a low dimensional basis)
# Pros: the data are much smoother and the measurement error is filtered
# Cons: I could have lost important information
basis.2 <- create.fourier.basis(rangeval=range(abscissa),nbasis=21)
data_W.fd.2 <- Data2fd(y = data_W,argvals = time,basisobj = basis.2)
plot.fd(data_W.fd.2)

# Choice 3: compromise between 1 and 2
basis.3 <- create.fourier.basis(rangeval=range(abscissa),nbasis=109)
data_W.fd.3 <- Data2fd(y = data_W,argvals = time,basisobj = basis.3)
plot.fd(data_W.fd.3)



# COEFFICIENTS
data_W.fd.1$coefs[howmanycoefs,whichobservations]



##### PERFORMING FPCA
###-------------------------------
# interpolated data (Choice 1)
plot.fd(data_W.fd.1,ylab='temperature')

pca_W.1 <- pca.fd(data_W.fd.1,nharm=5,centerfns=TRUE)
pca_W.1$varprop

# scree plot
# pca.fd computes all the 365 eigenvalues, but only the first 
# N-1=34 are non-null
plot(pca_W.1$values[1:35],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_W.1$values)[1:35]/sum(pca_W.1$values),xlab='j',ylab='CPV',ylim=c(0.8,1))
# how much explained by the first 3
cumsum(pca_W.1$values)[3]/sum(pca_W.1$values)

# plot eigenfunctions
par(mfrow = c(1,2))
plot(pca_W.1$harmonics[1,],col=1,ylab='FPC1',ylim=c(-0.1,0.08))
abline(h=0,lty=2)
plot(pca_W.1$harmonics[2,],col=2,ylab='FPC2',ylim=c(-0.1,0.08))

# plot of the FPCs as perturbation of the mean
# media <- mean.fd(data_W.fd.1)
# 
# plot(media,lwd=2,ylim=c(-25,20),ylab='temperature',main='FPC1')
# lines(media+pca_W.1$harmonics[1,]*sqrt(pca_W.1$values[1]), col=2)
# lines(media-pca_W.1$harmonics[1,]*sqrt(pca_W.1$values[1]), col=3)
# # variation in amplitude (more in winter than in summer)
# 
# plot(media,lwd=2,ylim=c(-20,20),ylab='temperature',main='FPC2')
# lines(media+pca_W.1$harmonics[2,]*sqrt(pca_W.1$values[2]), col=2)
# lines(media-pca_W.1$harmonics[2,]*sqrt(pca_W.1$values[2]), col=3)
# # temperate climate or not

# Command of the library fda that automatically does these plots
howmany <- c(1,2)
par(mfrow=howmany) # c(1,2) 
plot(pca_W.1, nx=100, pointplot=TRUE, harm=howmany, expand=0, cycle=FALSE)


# scatter plot of the scores
par(mfrow=c(1,2))
plot(pca_W.1$scores[,1],pca_W.1$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
points(pca_W.1$scores[35,1],pca_W.1$scores[35,2],col=2, lwd=4) # higlightin a particular object (35)

# observe single observation against and the mean
media <- mean.fd(data_W.fd.1)
day1 <- t(as.matrix(data[1,]))
plot.fd(data_W.fd.1[1,])
points(abscissa,day1) # time = abscissa
lines(media,lwd=2)


# plot(pca_W.1$scores[,1],pca_W.1$scores[,2],type="n",xlab="Scores FPC1",
#      ylab="Scores FPC2",xlim=c(-400,250))
# text(pca_W.1$scores[,1],pca_W.1$scores[,2],dimnames(data_W)[[2]], cex=1)

