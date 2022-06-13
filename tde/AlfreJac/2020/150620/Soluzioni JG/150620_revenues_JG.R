revenues <- read.table("revenues.txt", header = T)

library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics

coordinates(revenues) <- c('x','y')

x11()
bubble(revenues,'revenue',do.log=FALSE,key.space='bottom')
x11()
hist(revenues$revenue, breaks=16, col="grey", main='Histogram of rev', prob = TRUE, xlab = 'Zn')
x11()
xyplot(revenues$revenue ~ revenues$distance, as.data.frame(revenues))
#we have clearly negative correlation


svgm <- variogram(revenue ~ population, revenues)

x11()
plot(svgm, main = 'Sample Variogram',pch=19)
#we check isotropy
x11()
plot(variogram(revenue ~ population, revenues, alpha = c(0, 45, 90, 135)),pch=19)
#we can assume isotropy


#we try different fits
v.fit1 <- fit.variogram(svgm, vgm(500, "Sph", 2000))
x11()
plot(svgm, v.fit1, pch = 3)
v.fit1

v.fit2 <- fit.variogram(svgm, vgm(500, "Sph", 2000, 15))
x11()
plot(svgm, v.fit2, pch = 3)
v.fit2
#we do not have nugget

g <- gstat(formula = revenue ~ population, data = revenues, model = v.fit1)
g

new_a0 <- as.data.frame(cbind(x = 0, y = 0,population = 0))
coordinates(new_a0) <- c('x','y')

a0 <- predict(g, new_a0, BLUE = TRUE)$var1.pred# -> a0 = -29.02746

new_a1 <- as.data.frame(cbind(x = 0, y = 0,population = 1))
coordinates(new_a1) <- c('x','y')
a1 <- predict(g, new_a1, BLUE = TRUE)$var1.pred - a0# -> a1 = 0.0209

mod <- lm(population ~ distance, data = revenues)
summary(mod)
shapiro.test(mod$residuals)
x11()
plot(mod$fitted.values, mod$residuals/summary(mod)$sigma)
#so we assumed
#population = beta0 + beta1*distance + eps,  eps~ N(0,sigma^2)

s_duomo <- as.matrix(c(x = 514711.6, y = 5033903.0))
s_0 <- as.matrix(c(514703.8, 5035569.3))
d_0 <- norm(s_0 - s_duomo)
d_0bis <- sqrt((s_0[1]-s_duomo[1])^2 + (s_0[2]-s_duomo[2])^2)
p_0 <- predict(mod, data.frame(distance = d_0bis))
p_0bis <- mod$coefficients[1] + mod$coefficients[2]*d_0

z_new <- data.frame(cbind(x = s_0[1], y = s_0[2], population = p_0))
coordinates(z_new) <- c('x','y')
r_0 <- predict(g, z_new)$var1.pred# -> a1 = 0.0209

var <- predict(g, z_new)$var1.var
#no, this variance is defined by exploiting the covariance structure, estimated by fitting the sample variogram, but it doesn't account for the uncertainty over this estimation itself