traffic <- read.table("traffic.txt", header = T)
data <- t(traffic)
library(fda)
m <- 4           # spline order 
degree <- m-1    # spline degree 

abscissa <- 1:24
basis.1 <- create.bspline.basis(rangeval=c(1,24),nbasis=15,norder = m)

x11()
matplot(data,type='l',main='traffic',xlab='hour',ylab='volume')

data.fd <- Data2fd(y = data,argvals = abscissa, basisobj = basis.1)

x11()
plot.fd(data.fd, main = 'Smooted curves',xlab='hour',ylab='volume')

coeffs_day1 <- data.fd$coefs[1:3,1]

#PCA
pc <- pca.fd(data.fd,nharm=3,centerfns=TRUE)
pc
prop <- pc$values[1:3]/sum(pc$values)
cum_prop <- cumsum(prop)

x11()
layout(cbind(1,2))
plot(pc$values/sum(pc$values), pch = 5)

plot(cumsum(pc$values)/sum(pc$values), type = 'l', ylim = c(0,1))
abline(h = 0.8, col = 'blue')


x11()
layout(rbind(1,2,3))
plot(pc$harmonics[1,],col=1,ylab='FPC1')
plot(pc$harmonics[2,],col=2,ylab='FPC2')
plot(pc$harmonics[3,],col=2,ylab='FPC3')


#In order to interpret the results we consider:
media <- mean.fd(data.fd)

x11()
plot(media,lwd=2,ylab='volume',main='FPC1')
lines(media+pc$harmonics[1,]*sqrt(pc$values[1]), col=2)
lines(media-pc$harmonics[1,]*sqrt(pc$values[1]), col=3)
# Variation in flow, higher during the day and slightly lower during the night, wrt to the average, for positive scores

x11()
plot(media,lwd=2,ylab='volume',main='FPC2')
lines(media+pc$harmonics[2,]*sqrt(pc$values[2]), col=2)
lines(media-pc$harmonics[2,]*sqrt(pc$values[2]), col=3)
# days less trafficked in the morning wrt the average, but with more traffic from 10am -> 5 am

x11()
plot(media,lwd=2,ylab='volume',main='FPC3')
lines(media+pc$harmonics[3,]*sqrt(pc$values[3]), col=2)
lines(media-pc$harmonics[3,]*sqrt(pc$values[3]), col=3)
#more used during the afternoon wrt the average


#we can just limit ourself to the first principal component, which explains 90% of the variability by itself
#It might seem that the phenomenon can be just described in terms on when this road is used more than the average during each day

x11()
plot(1:30, pc$scores[,1])
#so during working days this road is more trafficked during day-time, while on weekend it's more used during night-times