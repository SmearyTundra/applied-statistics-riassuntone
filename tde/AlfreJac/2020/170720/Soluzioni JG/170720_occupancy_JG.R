occupancy <- read.table("occupancy.txt", header = T)
library(MASS)
attach(occupancy)

data <- occupancy[,1:2]

load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")


#we first of all need to check gaussianity
mcshapiro.test(data[which(X == 0),])
mcshapiro.test(data[which(X == 1),])
#the two groups are gaussian, we surely can use QDA, can we use LDA?
#covariance structure
S0 <- cov(data[which(X == 0),])
S1 <- cov(data[which(X == 1),])
S0
S1
var.test(data[which(X == 0),1],data[which(X == 1),1])
var.test(data[which(X == 0),2],data[which(X == 1),2])
#we cannot assume the same covariance structure, we must proceed with QDA

p = c(15,9)/24
q <- qda(data, X, prior = p)

h <- seq(min(Humidity), max(Humidity), length.out = 200)
c <- seq(min(CO2), max(CO2), length.out = 200)
hc <- expand.grid(Humidity = h, CO2 = c)
z <- predict(q,hc)$post
z1 <- z[,1] - z[,2]

x11()
plot(data)
points(data[which(X == 0),], col = 'red')
points(data[which(X == 1),], col = 'blue')
contour(h, c, matrix(z1, 200), levels=0, drawlabels=F, add=T)


#in order to compute a miningful APER we need to adjust it with the prior probabilities
Qda.m <- predict(q, data)
table(class.true=X, class.assigned=Qda.m$class)

APER  <- (1*p[1]/40+5*p[2]/60)
APER

z_new <- data.frame(Humidity = 26, CO2 = 9)
predict(q, z_new)$post
#it predicts it as an occupied room


library(class)
k <- knn(train = data, test = hc, cl = X, k = 5)

knn.class <- (k == 0)+0 
levels(as.factor(knn.class))

x11()
plot(data)
points(data[which(X == 0),], col = 'red')
points(data[which(X == 1),], col = 'blue')
contour(h, c, matrix(knn.class, 200), levels=0.5, drawlabels=F, add=T)
contour(h, c, matrix(z1, 200), levels=0, drawlabels=F, add=T, col = 'gray')



k$APER
predict(k, data)
kerr <- knn(train = data, test = data, cl = X, k = 5)
errs <- sum(as.numeric(kerr != X))
APERknn <- errs/dim(data)[1]
APERknn


errs_corr <- sum(ifelse(kerr != X,ifelse(X==0, p[1]/length(which(X==0)), p[2]/length(which(X==1))),0))

#the QDA has a smaller error rate, which moreover takes into account the prior probabilities
#even if we adjust the error estimate wrt prior probabilities, QDA is best