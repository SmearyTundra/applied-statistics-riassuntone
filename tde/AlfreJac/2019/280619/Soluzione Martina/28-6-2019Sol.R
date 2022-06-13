# EXAM 28-6-2019

#######################################
############# Exercise 1 ##############
#######################################

terrassa <- read.table('terrassa.txt', header=T)
girona <- read.table('girona.txt', header=T)

# a) Perform a statistical test of level 95% to verify if the mean evaluations in the two cities 
# differ. State and verify the model assumptions.

# The distribution of two different tapas are independent -> two differernt gaussian population!
# I can't do the dataset of differences

t1 <- terrassa
t2 <- girona

n1 <- dim(t1)[1] # n1=35
n2 <- dim(t2)[1] # n2=35
p  <- dim(t1)[2] # p=2


# we compute the sample mean, covariance matrices and the matrix Spooled

t1.mean <- sapply(t1,mean)
t2.mean <- sapply(t2,mean)
t1.cov  <-  cov(t1)
t2.cov  <-  cov(t2)
Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2)
# we compare the matrices
list(S1=t1.cov, S2=t2.cov, Spooled=Sp)
par(mfrow=c(1,2))
image(t1.cov, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(t1.cov,t2.cov), (0:100)/100, na.rm=TRUE))
image(t2.cov, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(t1.cov,t2.cov), (0:100)/100, na.rm=TRUE))
# same covariance structure!


# Test H0: mu1 == mu2  vs  H1: mu1 != mu2
# i.e.,
# Test H0: mu1-mu2 == c(0,0)  vs  H1: mu1-mu2 != c(0,0)

alpha   <- .05
delta.0 <- c(0,0)
Spinv   <- solve(Sp)

T2 <- n1*n2/(n1+n2) * (t1.mean-t2.mean-delta.0) %*% Spinv %*% (t1.mean-t2.mean-delta.0)

cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
T2 < cfr.fisher # TRUE: no statistical evidence to reject H0 at level 1%

P <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
P # 0  

# I can reject H0! The mean evaluation in the two cities differ.

# b) Interpret the results of the test at point (a) through two Bonferroni intervals of global 
# level 95% for appropriate differences in the mean. Comment the result

k <- p  # 2
cfr.t <- qt(1-alpha/(2*k),n1+n2-2)

IC.T2.1 <- c(t2.mean[1]-t1.mean[1]-sqrt(cfr.t*Sp[1,1]*(1/n1+1/n2)), t2.mean[1]-t1.mean[1]+sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)) )
IC.T2.2 <- c(t2.mean[2]-t1.mean[2]-sqrt(cfr.t*Sp[2,2]*(1/n1+1/n2)), t2.mean[2]-t1.mean[2]+sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)) )
IC.T2 <- rbind(IC.T2.X1, IC.T2.X2)
dimnames(IC.T2)[[2]] <- c('inf','sup')                        
IC.T2


# no interval contains the zero value, so I can confirm that on average there is a difference 
# between the city ratings

# c) Is there statistical evidence to state that, at level 95%, the average∗ evaluations of Girona’s 
# tapas are in mean higher than those of Terrassa’s tapas?

# Dataset of the averages

av.terrassa <- (terrassa[,1] + terrassa[,2])/2
av.girona <- (girona[,1] + girona[,2])/2

# My test: H0: mu(av.girona) >= mu(av.terrassa) vs H1: H0^C
#          H0: mu(av.girona) - mu(av.terrassa) >= 0 vs H1: H0^C

t1 <- av.terrassa
t2 <- av.girona

# Gaussianity assumption:
shapiro.test(t1)
shapiro.test(t2)
# ok

t1.cov  <-  var(t1)
t2.cov  <-  var(t2)
Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2)
# we compare the matrices
list(S1=t1.cov, S2=t2.cov, Spooled=Sp)
var.test(t1,t2) #ok!

t.test(t2,t1, alternative = 'less', var.equal = TRUE)

# pvalue of the test is 1, so I accept H0 and I can assess that there is statistical
# evidence to say that the evaluations of girona is higher than the one in terrassa.


#######################################
############# Exercise 2 ##############
#######################################

prof <- read.table('profiling.txt', header=T)

# a) Build a classifier for the variable type of user based on the available quantitative features.
# Report the model for the data, the estimates of its parameters (means and covariances), 
# the priors within the groups and verify the model assumptions. Report a qualitative plot of 
# the classification regions.

attach(prof)

# Look at the data:
plot(prof[,1], prof[,2], pch=19, col=c('blue','red')[factor(prof[,3])], xlab='t1', ylab='t2')

t <- prof[,1:2]
# Verify the assumptions for LDA:
# 1) normality (multivariate) within the groups
mcshapiro.test(t[which(prof$type=='tourist'),]) # pvalue=0.7
mcshapiro.test(t[which(prof$type=='resident'),]) # pvalue=0.13

# 2) equal variances
S1 <- cov(t[which(prof$type=='tourist'),])
S2 <- cov(t[which(prof$type=='resident'),])
par(mfrow=c(1,2))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
S1
S2

# I cannot use LDA approach
# I use qda

prof.lda <- qda(t,prof$type)
prof.lda

# Priors:
prof.lda$prior
# resident   tourist 
# 0.7132146 0.2867854 

# Plot of the classifcation regions
Lda.prof <- predict(prof.lda, t)

plot(t, main='t', xlab='t1', ylab='t2', pch=20)
points(t[which(prof$type=='tourist'),], col='red', pch=20)
points(t[which(prof$type=='resident'),], col='green', pch=20)
legend("topright", legend=levels(factor(prof$type)), fill=c('red','green'), cex=.7)

points(prof.lda$means, pch=4,col=c('red','green') , lwd=2, cex=1.5)

x  <- seq(min(prof[,1]), max(prof[,1]), length=200)
y  <- seq(min(prof[,2]), max(prof[,2]), length=200)
xy <- expand.grid(t1=x, t2=y)

z  <- predict(prof.lda, xy)$post
z1 <- z[,1] - pmax(z[,2])  
z2 <- z[,2] - pmax(z[,1])     

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

# b) Compute the APER of the classifier.

errors <- (Lda.prof$class != prof$type)
sum(errors)
length(prof$type)

APER   <- sum(errors)/length(prof$type)
APER # 0.0562

# c) How would you profile a new user with t1 = 35min and t2 = 3min?
x <- data.frame(t1 = 35, t2 = 3)

predict(prof.lda, x)$posterior

# I can assess that this profile is a tourist


#######################################
############# Exercise 3 ##############
#######################################

airport <- read.table('airport.txt', header=T)
attach(airport)

# a) Estimate the parameters of the model ({β0,g,β1,g,σ}). Verify the model assumptions.

# Dummies:
d1 <- ifelse(time.of.the.day=='6-10',1,0) # time 6-10
d2 <- ifelse(time.of.the.day=='11-15',1,0) # time 11-15
# If either = 0 -> time 16-20

model <- lm(duration ~ d1 + d2 + distance:d1 + distance:d2 + distance)
summary(model)

# Coefficients:  (group1 = 6-10, group2 = 11-15, group3 = 16-20)
coeffs <- data.frame(beta0.1 = model$coefficients[1] + model$coefficients[2],
                     beta0.2 = model$coefficients[1] + model$coefficients[3],
                     beta0.3 = model$coefficients[1],
                     beta1.1 = model$coefficients[4] + model$coefficients[5],
                     beta1.2 = model$coefficients[4] + model$coefficients[6],
                     beta1.3 = model$coefficients[4],
                     sigma = sqrt(sum(residuals(model)^2)/model$df))
coeffs

# Verify the model assumptions:
# Gaussianity on the residuals:
shapiro.test(model$residuals) # pvalue = 0.08 -> ok
# Homoschedasticity:
plot(model) # no particular pattern -> ok

# b) Perform two statistical tests – each at level 1% – to verify if
# - there is a significant dependence of the mean duration on the time of the day:
linearHypothesis(model, rbind(c(1,0,0,0,0,0),
                              c(0,1,0,0,0,0),
                              c(0,0,1,0,0,0)), c(0,0,0))
# there is statistical evidence to assess that the 3 factors are not significant for my model!
# I can remove these
model1 <- lm(duration ~ distance:d1 + distance:d2 + distance)
summary(model1)

# - there is a significant dependence of the mean duration on the distance traveled.
linearHypothesis(model1, rbind(c(0,1,0,0),
                              c(0,0,1,0),
                              c(0,0,0,1)), c(0,0,0))
# at least one is significant:
# looking at one at the time tests in the summary:
summary(model1)
# I can remove distance:d2:
model2 <- lm(duration ~ distance:d1 + distance)
summary(model2)

# OK! all the other terms are significant

# c) Based on tests (b) or any other test deemed relevant, reduce the model and update the 
# model parameters.

model.reduce <- lm(duration ~ distance:d1 + distance)
summary(model.reduce)

coeffs.new <- data.frame(beta0 = model$coefficients[1], 
                     beta1.1 = model$coefficients[2] + model$coefficients[3],
                     beta1.3 = model$coefficients[2],
                     sigma = sqrt(sum(residuals(model.reduce)^2)/model.reduce$df))
coeffs.new

# d) You have a flight from Barcelona airport at 10:30 a.m., and you want to be at the airport at 
# least 1 hour before the flight departure. At the bus station in front of your hotel in Terrassa 
# (distance of 57 km from the airport), the bus is scheduled to depart every 30 minutes from 6 
# a.m. to 20:30 p.m.. What time would you take the bus to be on time with probability 99%?

Z0.new <- data.frame(distance=57, d1=1) 

# Pred. int. for the new obs
Pred <- predict(model.reduce, Z0.new, interval='prediction', level=1-0.01)  
Pred

# fit      lwr      upr
# 112.4163 97.77151 127.0611

# the duration of the trip is ~ 112 mins, so if I want to be at the airport at 9.30, I must take 
# the bus at 7.30 a.m. from Terrassa (to stay calm I would take the one at 7 :P)


#######################################
############# Exercise 4 ##############
#######################################

montserrat <- read.table('montserrat.txt', header=T)

v.f <- function(x, ...){100-cov.spatial(x, ...)} # for graphical pourposes
v.f.est<-function(x,C0, ...){C0-cov.spatial(x, ...)} # for graphical pourposes

coordinates(montserrat) <- c('x','y')
bubble(montserrat,'speed',do.log=TRUE,key.space='bottom')

# a) First variogram:
svgm <- variogram(speed ~ 1, montserrat)
plot(svgm, main = 'Sample Variogram',pch=19)

vgm1 <- variogram(speed ~ distance, montserrat)
plot(vgm1, main = 'Second variogram', pch=19)

# The second variogram seems to increase linearly and than stabilized -> is more in line with
# the stationarity assumption, so I choose it!

# b) Fit to the empirical variogram chosen at point (a), a spherical model without nugget, via 
# weighted least squares. Report the estimates of sill and range. Comment the results.

v.fit1 <- fit.variogram(vgm1, vgm(1, "Sph", 10))
v.fit1

# estimate of sill and range:

#     sill    range
# 8.052468 28.32964

# c) Estimate, via Generalized Least Squares, the parameter(s) a of the model chosen at point (a)
z0 <- montserrat[1,]
z0$distance <- 0
g.t <- gstat(formula = speed ~ distance, data = montserrat, model = v.fit1)
# Estimated a0

a0 <- predict(g.t, z0, BLUE = TRUE)$var1.pred
a0 # 49.15834

# estimated a1
z0$distance <- 1
a1 <- predict(g.t, z0, BLUE = TRUE)$var1.pred - a0
a1 # -0.1019893

# d) Predict the wind speed at the top of the mountain, s0 = (402476, 4605558). Report the 
# associated prediction variance.
s0.new <- as.data.frame(matrix(c(402476,4605558,0),1,3))
names(s0.new) <- c('x','y','distance')
coordinates(s0.new) <- c('x','y')
predict(g.tr, s0.new)

# coordinates       var1.pred  var1.var
# (402476, 4605558)  52.49582  2.128486










