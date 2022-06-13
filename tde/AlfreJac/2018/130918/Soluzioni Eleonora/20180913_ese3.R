# EXAM 13/09/2018 ex 3

setwd("D:/Dati/Ele/poli/anno4/applied statistics/OldExams/tde")
load("D:/Dati/Ele/poli/anno4/applied statistics/R lab/LAB_5/mcshapiro.test.RData")

library(MASS)

# The file Sailing.txt collects the data on the daily values of consumed water
# [l/day] and sailing time [min/day] for 120 sailing cruises in Croatia, in
# August 2018. It also reports whether the sailboat was occupied by expert
# sailors (seadog) or by inexperienced vacationers (vacationer ). It has been
# estimated that in August, on average, only 20% of sailboats are occupied by
# expert sailors.

sailing <- read.table('Sailing.txt', header=T)
head(sailing)
attach(sailing)

# a) Based on the available features, build two Bayes classifiers, A and B, for
#    the kind of sailboat' occupation (vacationer, seadog), by assuming that:
#    A. the two populations are Gaussian with the same covariance structure;
#    B. the two populations are Gaussian with different covariance structures.
#    For each classifier report a qualitative plot of the classification regions,
#    and the estimated posterior probability associated with the first
#    observation (water = 32.08, sailing.time = 82.69 ).

# Prior probabilities
p1 <- 0.2
p2 <- 0.8

group.name <- factor(type, labels = c('seadog', 'vacationer'))
sail <- sailing[,1:2]
i1 <- which(type=='seadog')
i2 <- which(type=='vacationer')

# CLASSIFIER A:
# The two populations are Gaussian with the same covariance structure
# -> LDA

lda.A <- lda(sail, group.name, prior=c(p1,p2))

x11()
plot(sail, main='Sailing type A', xlab='water', ylab='sailing.time', pch=20)
points(sail[i1,], col='red', pch=20)
points(sail[i2,], col='green', pch=20)
legend("topright", legend=levels(group.name), fill=c('red','green'), cex=.7)

points(lda.A$means, pch=4,col=c('red','green') , lwd=2, cex=1.5)

x  <- seq(min(sailing[,1]), max(sailing[,1]), length=200)
y  <- seq(min(sailing[,2]), max(sailing[,2]), length=200)
xy <- expand.grid(water=x, sailing.time=y)

z  <- predict(lda.A, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - z[,2]  # P_1*f_1(x,y)-P_2*f_2(x,y)  
z2 <- z[,2] - z[,1]  # P_2*f_2(x,y)-P_1*f_1(x,y) 

# Plot the contour line of level (levels=0) of z1, z2: 
# P_i*f_i(x,y)-P_j*f_j(x,y)=0 i.e., boundary between R.i and R.j 
# where j realizes the max.
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

# Posterior probabilities
LDA.A <- predict(lda.A, sail)
LDA.A$posterior[1,]
# seadog     vacationer 
# 0.95565267 0.04434733
LDA.A$class[1] # seadog

# CLASSIFIER B:
# The two populations are Gaussian with different covariance structures
# -> QDA

qda.B <- qda(sail, group.name, prior=c(p1,p2))

x11()
plot(sail, main='Sailing type B', xlab='water', ylab='sailing.time', pch=20)
points(sail[i1,], col='red', pch=20)
points(sail[i2,], col='green', pch=20)
legend("topright", legend=levels(group.name), fill=c('red','green'), cex=.7)

points(qda.B$means, pch=4,col=c('red','green') , lwd=2, cex=1.5)

x  <- seq(min(sailing[,1]), max(sailing[,1]), length=200)
y  <- seq(min(sailing[,2]), max(sailing[,2]), length=200)
xy <- expand.grid(water=x, sailing.time=y)

z  <- predict(qda.B, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - z[,2]  # P_1*f_1(x,y)-P_2*f_2(x,y)  
z2 <- z[,2] - z[,1]  # P_2*f_2(x,y)-P_1*f_1(x,y) 

# Plot the contour line of level (levels=0) of z1, z2: 
# P_i*f_i(x,y)-P_j*f_j(x,y)=0 i.e., boundary between R.i and R.j 
# where j realizes the max.
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

# Posterior probabilities
QDA.B <- predict(qda.B, sail)
QDA.B$posterior[1,]
# seadog vacationer 
# 0.96474706 0.03525294 
QDA.B$class[1] # seadog

# ___________________________________________________________________________

# b) Evaluate the performances of the classifiers A and B and identify the best
#    one.

prior <- c(p1,p2)
G <- 2

# CLASSIFIER A:
misc.A <- table(class.true=group.name, class.assigned=LDA.A$class)
APER.A <- 0
for (g in 1:G)
  APER.A <- APER.A + sum(misc.A[g,-g])/sum(misc.A[g,])*prior[g]
APER.A #0.06333333

# CLASSIFIER B:
misc.B <- table(class.true=group.name, class.assigned=QDA.B$class)
APER.B <- 0
for (g in 1:G)
  APER.B <- APER.B + sum(misc.B[g,-g])/sum(misc.B[g,])*prior[g]
APER.B #0.06

# The apparent error rate of  is less than the one of A, so the best classifier
# is B

# ____________________________________________________________________________

# c) How would you classify the occupants of a sailboat with daily consumed
#    water 35 l and daily sailing time 168 min?

new.data <- data.frame(water=35, sailing.time=168)

# CLASSIFIER A:
pred.A <- predict(lda.A, new.data)
pred.A$posterior
# seadog     vacationer
# 0.8061144  0.1938856
pred.A$class # seadog

# CLASSIFIER B:
pred.B <- predict(qda.B, new.data)
pred.B$posterior
# seadog     vacationer
# 0.1577169  0.8422831
pred.B$class # vacationer

points(35,168, pch=20, col='blue')
