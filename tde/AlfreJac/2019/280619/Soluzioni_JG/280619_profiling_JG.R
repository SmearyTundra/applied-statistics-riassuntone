profiling <- read.table("profiling.txt", header = T)
n <- dim(profiling)[1]
attach(profiling)
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
tour <- which(type=="tourist")
mcshapiro.test(profiling[tour, 1:2])
mcshapiro.test(profiling[-tour, 1:2])
#we can assume gaussianity in both populations
St <- cov(profiling[tour, 1:2])
Sr <- cov(profiling[-tour, 1:2])
St
Sr
#we cannot assume homoschedasticity, we proceed with QDA
#For which Rt = {x e Rp: log(pt) - 1/2log(det(SIGMAt) - 1/2(x - mut)(SIGMAt)^-1(x - mut) >= log(pr) - 1/2log(det(SIGMAr) - 1/2(x - mur)(SIGMAr)^-1(x - mur)}
library(MASS)
#So 
q <- qda(type ~ t1 + t2)
q

x <- seq(min(t1), max(t1), length = 200)
y <- seq(min(t2), max(t2), length = 200)
xy <- expand.grid(t1 = x, t2 = y)
z <- predict(q,xy)$post
z1 <- z[,1]-z[,2]

x11()
plot(profiling[,1:2], main = "Induced partition - QDA")
points(profiling[tour, 1:2], col = 'red')
points(profiling[-tour, 1:2], col = 'blue')
legend("topright", legend=levels(as.factor(type)), fill=c('blue','red'), cex=.7)
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  



#APER
q.APER <- predict(q, profiling[,1:2])$class
t <- table(true = type, assigned = q.APER)
t
APER <- (t[1,2]+t[2,1])/n
APER

z_new <- data.frame(t1 = 35, t2 = 3)
predict(q, z_new)$post
#almost surely a tourist