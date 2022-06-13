sansiro <- read.table("sansiro.txt",header = T)
attach(sansiro)
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
mcshapiro.test(sansiro[which(opera == "regular"),1:2])
mcshapiro.test(sansiro[which(opera == "suspend"),1:2])
mcshapiro.test(sansiro[which(opera == "cancel"),1:2])
#we can assume gaussianity between groups
S_r <- cov(sansiro[which(opera == "regular"),1:2])
S_s <- cov(sansiro[which(opera == "suspend"),1:2])
S_c <- cov(sansiro[which(opera == "cancel"),1:2])
S_r
S_s
S_c
bartlett.test(precipitation ~ opera)
bartlett.test(temperature ~ opera)
#we can assume homoschedasticity

#We can hence use LDA
pr <- c(3, 100-8-3, 8)/100

library(MASS)
l <- lda(opera ~ temperature + precipitation, prior = pr)
l
#we can get here the group means
t <- seq(min(temperature), max(temperature), length = 200)
p <- seq(min(precipitation), max(precipitation), length = 200)
tp <- expand.grid(temperature = t, precipitation = p)

z <- predict(l, tp)$post

z1 <- z[,1] - pmax(z[,2], z[,3])
z2 <- z[,2] - pmax(z[,1], z[,3])

x11()
plot(sansiro[,1:2], main = "Partition LDA")
points(sansiro[which(opera == "regular"),1:2], col = 'green')
points(sansiro[which(opera == "cancel"),1:2], col = 'red')
points(sansiro[which(opera == "suspend"),1:2], col = 'blue')
contour(t, p, matrix(z1, 200), levels=0, drawlabels=F, add=T)
contour(t, p, matrix(z2, 200), levels=0, drawlabels=F, add=T)

#computation of APER
pred.class <- predict(l, sansiro[,1:2])$class
table(class.true = opera, class.assigned = pred.class)
APER <- pr[1]*1/15 + pr[2]*0 + pr[3]*5/15
APER

z_new <- data.frame(temperature = 28, precipitation = 12)
predict(l, z_new)$posterior
#it will probably be suspended


#we can build a linear classifier which takes into account only precipitation and compare the APER
l_red <- lda(opera ~ precipitation, prior = pr)
l_red

z_red <- predict(l_red, tp)$post

z1_red <- z_red[,1] - pmax(z_red[,2], z_red[,3])
z2_red <- z_red[,2] - pmax(z_red[,1], z_red[,3])

x11()
plot(sansiro[,1:2], main = "Partition LDA")
points(sansiro[which(opera == "regular"),1:2], col = 'green')
points(sansiro[which(opera == "cancel"),1:2], col = 'red')
points(sansiro[which(opera == "suspend"),1:2], col = 'blue')
contour(t, p, matrix(z1_red, 200), levels=0, drawlabels=F, add=T)
contour(t, p, matrix(z2_red, 200), levels=0, drawlabels=F, add=T)
contour(t, p, matrix(z1, 200), levels=0, drawlabels=F, add=T, col = 'grey')
contour(t, p, matrix(z2, 200), levels=0, drawlabels=F, add=T, col = 'grey')


pred.class_red <- predict(l_red, sansiro[,1:2])$class
table(class.true = opera, class.assigned = pred.class_red)
APER_red <- pr[1]*0/15 + pr[2]*0/15 + pr[3]*5/15
APER
APER_red
#even a better performance (on the training set)

#if we were to compare them via cross validation:
l.CV <- lda(opera ~ temperature + precipitation, prior = pr, CV = T)
table(class.true=opera, class.assignedCV=l.CV$class)
AER.CV <- pr[1]*1/15 + pr[2]*0/15 + pr[3]*5/15

l.CV_red <- lda(opera ~ precipitation, prior = pr, CV = T)
table(class.true=opera, class.assignedCV=l.CV_red$class)
AER.CV_red <- pr[1]*0/15 + pr[2]*0/15 + pr[3]*5/15

AER.CV
AER.CV_red


#we can neglect temperature