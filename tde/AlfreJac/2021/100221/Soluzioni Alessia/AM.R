######
# PROBLEMA 1
#####
library(car)

data <- read.table('shopping.txt',header = T)
head(data)
dim(data)
load("~/Desktop/IV anno/Applied statistics/mcshapiro.test.RData")
mcshapiro.test(data)
M <- sapply(data,mean)
S <- cov(data)
alpha <- 0.05
n <- dim(data)[1]
p <- dim(data)[2]
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)


# Center:
M
# Directions of the principal axes:
eigen(S/n)$vectors

# Length of the semi-axes of the ellipse:
r <- sqrt(cfr.fisher)
r*sqrt(eigen(S/n)$values) 

a <- c(-0.2,1,1)
#, c(0,1,0),c(0,0,1), c(1,1,1))
T2 <- cbind(inf = t(a)%*%M - sqrt(cfr.fisher*(t(a)%*%S%*%a)/n),
            center = t(a)%*%M , 
            sup = t(a)%*%M  + sqrt(cfr.fisher*(t(a)%*%S%*%a)/n))
T2

cfr.t<- qt(1-alpha,n-1)
T0 <- t(a)%*%M/sqrt(t(a)%*%S%*%a/n)
T0 > cfr.t
p <- 1-pt(T0, n-1)
p

######
# PROBLEMA 2
#####
data <- read.table('rice.txt', header=T)
dim(data)
head(data)
e <- dist(data, method='euclidean')
ec <- hclust(e, method='single') 


plot(ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
cluster.ec <- cutree(ec, k=3) # euclidean-complete:
cluster.ec
i1 <- which(cluster.ec== 1)
i2 <- which(cluster.ec== 2)
i3 <- which(cluster.ec== 3)

length(i1)
length(i2)
length(i3)

sapply(data[i1,], mean)
sapply(data[i2,], mean)
sapply(data[i3,], mean)

plot(data)
points(data[i1,], col= 'red')
points(data[i2,], col= 'blue')
points(data[i3,], col= 'green')
coph.ec <- cophenetic(ec)
ec <- cor(e, coph.ec)
ec

inf <- data[i3,]
load("~/Desktop/IV anno/Applied statistics/mcshapiro.test.RData")

mcshapiro.test(inf)
alpha <- 0.05
n <- dim(inf)[1]
k <- 2 # number of intervals I want to compute (set in advance)
cfr.t <- qt(1-alpha/(2*k),n-1)
x.mean <- sapply(inf,mean)
x.cov <- cov(inf)
Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf

CI.var <- cbind(inf = x.cov[1,1]*(n-1)/qchisq(1-alpha/6,n-1),
                              center = x.cov[1,1],
                sup = x.cov[1,1]*(n-1)/qchisq(alpha/6,n-1))
CI.var
CI.var <- cbind(inf = x.cov[2,2]*(n-1)/qchisq(1-alpha/6,n-1),
                center = x.cov[2,2],
                sup = x.cov[2,2]*(n-1)/qchisq(alpha/6,n-1))
CI.var
#####
#PROBLEMA 3
######
data <- read.table('landslides.txt',header = T)
head(data)
dim(data)
fm <- lm(rate~. , data)
summary(fm)
fm$coefficients

par(mfrow= c(2,2))
plot(fm)
shapiro.test(fm$residuals)

fm <- lm(rate~.-hardness , data)
summary(fm)
fm$coefficients

par(mfrow= c(2,2))
plot(fm)
shapiro.test(fm$residuals)

linearHypothesis(fm, c(0,0,1,-2), 0)
fm <- lm(rate~ rain+ I(2*coarse+fine) , data)
summary(fm)
fm$coefficients

par(mfrow= c(2,2))
plot(fm)
shapiro.test(fm$residuals)

z0 = data.frame(rain = 700,coarse=10,fine=8)
t(z0)%*%fm$coefficients
Conf <- predict(fm, z0, interval='confidence', level=1-0.01)  
Conf

######
# PROBLEMA 4 
######
data <- read.table('hotels.txt', header = T)
dim(data)
head(data)
library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics

## Functions for graphics 
v.f <- function(x, ...){100-cov.spatial(x, ...)}
v.f.est<-function(x,C0, ...){C0-cov.spatial(x, ...)}
coordinates(data) <- c('x','y')

hist(data$price, breaks=16, col="grey", main='Histogram of Zn', prob = TRUE, xlab = 'Zn')

bubble(data,'price',do.log=F, key.space='bottom')

v<- variogram(price ~ 1, data)
svgm
plot(svgm, main = 'Sample Variogram',pch=19)
plot(variogram(price ~ 1, data, alpha = c(0, 45, 90, 135)),pch=19)


v.fit <-vgm(600, "Sph", 1500, 1000)
v.fit <- fit.variogram(v, v.fit)
plot(v, v.fit, pch = 19)

s0.new=data[1,] # UTM coordinates
g.tr <- gstat(formula = price ~ 1, data = data, model = v.fit)
predict(g.tr, s0.new, BLUE = TRUE)
dummy <- rep(0,dim(data)[1])
dummy[which(data$winter== 'no' )] <- 1
data <- data.frame(data,dummy= dummy)
data <- data[,-6]
coordinates(data) <- c('x','y')
v<- variogram(price ~ 1 + dummy + distance + distance:dummy, data)
v
plot(v, main = 'Sample Variogram',pch=19)
plot(variogram(price ~ 1 + dummy + distance + distance:dummy, data, alpha = c(0, 45, 90, 135)),pch=19)


v.fit <-vgm(1000, "Sph", 1500, 250)
v.fit <- fit.variogram(v, v.fit)
plot(v, v.fit, pch = 19)

g.tr <- gstat(formula = price ~ 1 + dummy + distance + distance:dummy, data = data, model = v.fit)

Z0_for_params <- data[1,]
Z0_for_params$distance[1]=0
Z0_for_params$dummy = 0

predict(g.tr,Z0_for_params, BLUE = TRUE)
a01 <- predict(g.tr,Z0_for_params, BLUE = TRUE)$var1.pred

Z0_for_params$dummy = 1
predict(g.tr,Z0_for_params, BLUE = TRUE)
a02 <- a01 + predict(g.tr,Z0_for_params, BLUE = TRUE)$var1.pred

Z0_for_params$distance[1]=1
Z0_for_params$dummy = 0
predict(g.tr,Z0_for_params, BLUE = TRUE)
a11 <- predict(g.tr,Z0_for_params, BLUE = TRUE)$var1.pred - a01

Z0_for_params$distance[1]=1
Z0_for_params$dummy = 1
predict(g.tr,Z0_for_params, BLUE = TRUE)
a12 <- predict(g.tr,Z0_for_params, BLUE = TRUE)$var1.pred - a02


D.s0=sqrt( (342399.74 - 342362.58)^2 + (5072272.75-5072518.24)^2)
s0=as.data.frame(matrix(c(342399.74, 5072272.75,D.s0,0),1,4))
names(s0)=c('x','y','distance','dummy')
coordinates(s0)=c('x','y')

predict(g.tr, s0, BLUE = FALSE)
predict(g.tr, s0, BLUE = TRUE)
