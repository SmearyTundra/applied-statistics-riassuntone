rm(list = ls())
hotels = read.table("hotels.txt")
load("~/Desktop/Applied Stat - Lab/TdeIta/mcshapiro.test.RData")
library(car)
library(heplots)
library(Hotelling)
library(ICSNP)
library(rstatix)
library(MASS)
library(leaps)
library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics

# BLUE TRUE
#    - se vuoi stimare i coefficienti (GLS)
#    - not using completely the spacial information (se devi predirre qualcosa di indipendente)
# BLUE FALSE
#    - WLS
#    - se vuoi fare una predizione relativa a quel fenomeno

coordinates(hotels) <- c('x','y') # set the coordinates correctly

# bubble plot(obj,zcol,...)
# key.space=location of the key
bubble(hotels,'price',do.log=TRUE,key.space='bottom')



hist(hotels$price, breaks=16, col="grey", main='Histogram of price', prob = TRUE, xlab = 'Zn')
# highly skewed, transform to the log
# hist(log(hotels$price), breaks=16, col="grey", main='Histogram of log(Zn)', prob = TRUE, xlab = 'log(Zn)')

# scatterplot of log(zinc) with respect to distance from the river 
xyplot(price ~ distance, as.data.frame(hotels))

# how to assume isotropy
plot(variogram(log(sights) ~ 1 + log.chlorofill, data=data, alpha = c(0, 45, 90, 135)),pch=19)

v.t=variogram(price ~ winter + winter:distance , data=hotels)
plot(v.t,pch=19)





# list variogram models
vgm()
# "Exp" "Sph" are linear
#both spherical and exponential model have a linear behavior near the origin but exponential model has a faster growth than the spherical one


# sample variogram (binned estimator) + plot
v.no=variogram(price ~ 1 , data=hotels)
v.fit1 <- fit.variogram(v.no, vgm(4000, "Sph", 500,500))
		#vgm(sill, the value it approaches
		#	'type', Exp or Sph
		#	range, where is the elbow
		#	nugget jump at the origin, possibly omitted
		#	)
plot(v.no, v.fit1, pch = 3)
v.fit1

v.fit2 <- fit.variogram(v.t, vgm(1000, "Sph", 500,500))
plot(v.t, v.fit2, pch = 3)
v.fit2

# OK
g.no <- gstat(formula = price ~ 1 , data = hotels, model = v.fit1)
# UK
g.t <- gstat(formula = price ~ winter + winter:distance , data = hotels, model = v.fit2)

# predict hotels[1,]
# Estimate the mean: use the argument 'BLUE=TRUE' otherwise the observation
predict(g.no, hotels[1,], BLUE = TRUE)$var1.pred    

##### PREDICTION FOR GROUP 1  

s1 = hotels[1,]
w1 = predict(g.t, s1, BLUE = TRUE)$var1.pred
s2 = hotels[2,]
w2 = predict(g.t, s2, BLUE = TRUE)$var1.pred
a1 = (w2-w1)/(s2$distance - s1$distance)
a0 = w1 - a1*s1$distance
c(a0=a0,a1=a1)



datum=as.data.frame(t(matrix(c(342399.74,5072272.75,D.s0))))
names(datum)=c('x','y','distance')
coordinates(datum)=c('x','y')
datum = cbind(datum, winter='yes')
colnames(datum@data)[2] = "winter"

predict(g.t, datum, BLUE = TRUE)  ### BLUE TRUE perchè sto predicendo il prezzo del futuro

######### WE PASS TO ESTIMATE THE PARAM  (COEFF OF LINEAR MODEL)

a0 <- predict(g.tr, rv[1,], BLUE = TRUE)$var1.pred   #where variable = 0
a1 <- predict(g.tr, rv[6,], BLUE = TRUE)$var1.pred   #where variable = 1

b0 <- a0
b1 <- a1 - a0


###################
s1=hotels[1,]

w1= predict(g.t, s1,BLUE = TRUE)$var1.pred

s2=hotels[30,]

w2= predict(g.t, s2,BLUE = TRUE)$var1.pred

a1 = (w2-w1)/(s2$distance - s1$distance)
a0 = w1 - a1*s1$distance

param = c(a0,a1)
