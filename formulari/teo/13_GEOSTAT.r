# GEOSTATISTICS


library(gstat)
library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         

# Model assumptions

# - The ordinary Kriging model assumes II order stationariety ( E[z_s] = m for any s and Cov(z_s1,z_s2) = C(h) where h = norm(s_1 - s_2) ), and isotropy. We also assume that C(•) is known.

# - The universal Kriging model assumes z_s = m_s + delta_s for any s in D (domain) where m_s is called drift and describes the non constant spacial mean variation. Moreover we assume E[delta_s] = 0 for any s in D (so that E[z_s] = m_s) and that Cov(z_s1,z_s2) = Cov(delta_s1,delta_s2) for any pair. We also assume that C(•) is known and that m_s follows a linear model m_s(t) = sum_l=0^L a_l(t) f_l(s) for s in D and t in T, where f are known functions of s and a_l are coefficients independent from the spacial location.

# BLUE TRUE
#    - se vuoi stimare i coefficienti (GLS)
#    - not using completely the spacial information (se devi predirre qualcosa di indipendente)
# BLUE FALSE
#    - WLS
#    - se vuoi fare una predizione relativa a quel fenomeno

### READING COORDINATES
### -----------------------
coordinates(d) <- c('x','y')
bubble(d,'sights',do.log=TRUE,key.space='bottom')



### FITTING VARIOGRAM (ANISOTROPY AND NOT)
### -----------------------
v.gls <- variogram(log(sights) ~ 1+log.chlorofill,d)

# no anis
plot(variogram(log(sights) ~ 1+log.chlorofill,d),pch=19) 
# y anis
plot(variogram(log(sights) ~ 1+log.chlorofill, d, alpha = c(0, 45, 90, 135)),pch=19)

# doing fit
# VGM(PSILL, MODEL, RANGE, (NUGGET)) don't put nugget if you dont want it
#vgm(sill, the value it approaches
#	'type', Exp or Sph
#	range, where is the elbow
#	nugget jump at the origin, possibly omitted
#	)
v <- variogram(log(sights) ~ 1+log.chlorofill, d)
v.fit <- fit.variogram(v, vgm(0.5, "Exp", 50000))  
plot(v, v.fit, pch = 19)
v.fit


### CREATING GSTAT OBJECT
### -----------------------
g.no <- gstat(id = 'sights', formula = log(sights) ~ 1+log.chlorofill, data = d, model = v.fit)
g.no$model # MODEL ESTIMATED FOR DELTA(Si)



### ESTIMATING a0,a1 etc
### -----------------------
# Find smart combination, if stationary just the mean
x0 = data.frame(x=20,y=30,log.chlorofill=0)
coordinates(x0) <- c("x","y")
a0 <- predict(g.no, x0, BLUE = TRUE)$sights.pred #2.74
a0

x0 = data.frame(x=20,y=30,log.chlorofill=1)
coordinates(x0) <- c("x","y")
a1 <- predict(g.no, x0, BLUE = TRUE)$sights.pred - a0
a1

c(a0=a0,a1=a1)

# Smart combination for classical linear system
# y(si) = a0,g + a1,g · distance(si) + δ(si) g=1,2
# We have to estimate 4 coefficients.
s1 = hotels[1,]
w1 = predict(g.no, s1, BLUE = TRUE)$var1.pred
s2 = hotels[2,]
w2 = predict(g.no, s2, BLUE = TRUE)$var1.pred
a1 = (w2-w1) / (s2$distance - s1$distance)
a0 = w1 - a1 * s1$distance
c(a0=a0,a1=a1)


### PREDICTING
### -----------------------
x1 = data.frame(
	x=253844.8,
	y=385997.7,
	log.chlorofill=log.chlo)
coordinates(x1) <- c("x","y")
pr <- predict(g.no, x1)$sights.pred
pr

#C'è già sopra la var: no non è rappresentativa! Si usa un universal kriging quindi la varianza è ampiamente
#sottostimata, it does not account for the fact that sigma is unknown and it is estimated just running the algorithm (large underestimation of the uncertainty)