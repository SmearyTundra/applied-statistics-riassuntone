watertemp <- read.table("watertemp.txt", header = T)
library(fda)
data <- t(watertemp[,-length(watertemp)])#we eliminate the last column (about the zone) and transpose the dataset

x11()
matplot(data)

time <- 1:365
basis <- create.fourier.basis(rangeval=c(0,365),nbasis=45)

data.fd <- Data2fd(y = data,argvals = time,basisobj = basis)
x11()
plot.fd(data.fd)


coeff_point_a <- data.fd$coefs[1:3,1:2]

#PCA
pc <- pca.fd(data.fd, nharm = 5, centerfns = T)
var <- pc$varprop    # variability along the single PC
cumvar <- cumsum(pc$varprop)[1:5]  # cumulative variability

x11()
layout(rbind(1,2,3))
plot(pc$harmonics[1,],col=1,ylab='FPC1')
abline(h=0,lty=2)
plot(pc$harmonics[2,],col=2,ylab='FPC2')
plot(pc$harmonics[3,],col=2,ylab='FPC3')

x11()
plot(pc$values,xlab='j',ylab='Eigenvalues')
x11()
plot(cumvar,xlab='j',ylab='CPV',ylim=c(0.8,1))

#To interpret the PC it's useful to see the difference on 1 sd wrt the mean

media <- mean.fd(data.fd)

x11()
plot(media,lwd=2,ylim = c(5, 30), ylab='temperature',main='FPC1')
lines(media+pc$harmonics[1,]*sqrt(pc$values[1]), col=2)
lines(media-pc$harmonics[1,]*sqrt(pc$values[1]), col=3)

#PC1 -> Practically a variation in the mean temperature, still following the trend of cold winters and warmer summers
#the higher the hotter

x11()
plot(media,lwd=2,ylim = c(5, 30), ylab='temperature',main='FPC2')
lines(media+pc$harmonics[2,]*sqrt(pc$values[2]), col=2)
lines(media-pc$harmonics[2,]*sqrt(pc$values[2]), col=3)

#PC2 -> a time translation of the trend -> the higher the sooner the heat wave starts

x11()
plot(media,lwd=2,ylim = c(5, 30), ylab='temperature',main='FPC3')
lines(media+pc$harmonics[3,]*sqrt(pc$values[3]), col=2)
lines(media-pc$harmonics[3,]*sqrt(pc$values[3]), col=3)

#PC3 -> temperature excursions throughout the year -> the higher the smaller the variation between temperature in winter/summer is

x11()
plot(pc$scores[,1], pc$scores[,2])
points(pc$scores[,1][which(watertemp$Zone == "Surface")], pc$scores[,2][which(watertemp$Zone == "Surface")], col = 'red')
points(pc$scores[,1][which(watertemp$Zone == "Medium")], pc$scores[,2][which(watertemp$Zone == "Medium")], col = 'green')
points(pc$scores[,1][which(watertemp$Zone == "Deep")], pc$scores[,2][which(watertemp$Zone == "Deep")], col = 'blue')

#so deep waters tend to have lower temperatures with less variability on when the temperature change
#surface waters have a higher temperature and much more variability in terms of time-behaviour (summer heat pick may come much later/sooner)
#medium waters do not differ much from deep waters in terms of time-behaviour, but are significantly warmer


#We can reduce the dimensionality by just considering the first two PC, since they still account for 98% of the variability of our data