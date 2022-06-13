montserrat <- read.table("montserrat.txt", header = T)
library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics
coordinates(montserrat)<- c('x', 'y')

x11()
xyplot(speed ~ distance, as.data.frame(montserrat))

m1.svg<-variogram(speed ~ 1, montserrat)
m2.svg<-variogram(speed ~ distance, montserrat)
x11()
plot(m1.svg, main = "Sample Variogram M1",pch = 19)
x11()
plot(m2.svg, main = "Sample Variogram M2",pch = 19)
#the variogram of M1 seems to be divergent, while M2 seems to stabilize
#M2 is more robust with the hypothesis of stationarity, which entails a finite variance
v.fit <- fit.variogram(m2.svg, vgm(7, "Sph", 25))
v.fit

x11()
plot(m2.svg, v.fit, pch = 19)

g <- gstat(formula = speed ~ distance, data =montserrat, model = v.fit)
res.vgm <- variogram(g)
x11()
plot(res.vgm, main = "Residual Variogram M2",pch = 19)

res.fit <- fit.variogram(res.vgm, vgm(7, "Sph", 25))
res.fit
v.fit
#the assumptions of the model are verified

z_0 <- data.frame(x=0,y=0,distance=0)
coordinates(z_0) <- c('x','y')
a0 <- predict(g, z_0, BLUE = T)$var1.pred

z_1 <- data.frame(x=0,y=0,distance=1)
coordinates(z_1) <- c('x','y')

a1 <-  predict(g, z_1, BLUE = T)$var1.pred - a0

z_top <- data.frame(x = 402476, y = 4605558, distance = 0)
coordinates(z_top)<- c('x', 'y')
#we are looking for the actual prediction, not the trend component at that point
predict(g, z_top, BLUE = FALSE)

