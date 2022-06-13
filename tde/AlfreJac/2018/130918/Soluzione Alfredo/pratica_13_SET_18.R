setwd("~/Documents/Polimi/APPSTATS/")

library(mvtnorm)
library(mvnormtest)

iamg <- read.table("IAMG.txt")
head(iamg)
mcshapiro.test(iamg)
alpha <- 0.05
X <- sapply(iamg, mean)
X
S <- cov(iamg)
n <- dim(iamg)[1]
p <- dim(iamg)[2]
"Apply Hotelling's theorem, plausibility of X as mean value"
T2 <- n * t(X) %*% solve(S) %*% X
raggio <- (n-1)*p/(n-p) * qf(1-alpha, p, n-p)

# directions
eigen(S)[,1]
eigen(S)[,2]
eigen(S)[,3]

# lengths of semi-axes
sqrt(eigen(S)$values[1] * 1/n * raggio)
sqrt(eigen(S)$values[2] * 1/n * raggio)
sqrt(eigen(S)$values[3] * 1/n * raggio)


## b) simulatenous intervals
C <- matrix(c(1,0,0,
              0,1,0,
              0,0,1), byrow = T, nrow = 3, ncol=3)

C %*% X + sqrt(diag(C %*% S %*% t(C)) * 1/n) * sqrt(raggio)
C %*% X - sqrt(diag(C %*% S %*% t(C)) * 1/n) * sqrt(raggio)

####
### C)
proportions <- iamg$Talk / iamg$Registered
proportions
phat <- mean(proportions)
"Perform the transformation suggested from book from 
transforms to recover normality"
props.tr <- 1/2 *log(proportions/(1-proportions))
hist(props.tr)
shapiro.test(props.tr)
"So now we have samples from a normal distribution"
mu0 <- 1/2 *log(0.9/(1-.9))
"Perform a one sided tst
H0: mu = mu0 HA"
?t.test
"Perform test with H1 mu0 > 0.9, as we want to verify it is 'only' 0.9"
t.test(props.tr, mu = mu0, alternative="greater") 
"Reject the null hypothesis in favor of alternative that the mean
is less than 0.9"
# BONUS: perform the test generalizing from T2
phat2 <- mean(props.tr)
S <- var(props.tr)
T2.2 <- t(phat2-mu0) %*% S**-1 %*% (phat2-mu0)
cfr.F <- (n-1)*1/(n-1)*qf(1-alpha,1, n-1)
T2.2 < cfr.F
1- pf(T2.2, 1, n-1)  # p value more or less the same no 1-sided


###########
## PROBLEM 2
########|
wait <- read.table("Waiting.txt")
head(wait)


###
## a)
###
"ANOVA: y = mean + taui + beta.l + gamma. + eps.ij "
attach(wait)
course <- as.factor(course)
city <- as.factor(city)
levels(course)
levels(city)

bartlett.test(waiting, paste0(as.character(course),as.character(city)))
an <- aov(waiting ~ course + city + course:city)
summary(an)
shapiro.test(an$residuals)

## b)
an2 <- aov(waiting ~ course + city )
summary(an2)
# make a reduced model since neither interaction nor city are significant 
an3 <- aov(waiting ~ course )
summary(an3)
shapiro.test(an3$residuals)
bartlett.test(waiting, course)

##
## c)
##
k <- 3*(3-1)/2 + 1 # one for differences, one each within group variance  
DF <- an3$df.residual
N <- length(waiting)
Sp <- (t(an3$residuals) %*%an3$residuals) / DF
mediae <- tapply(waiting, course, mean)
an3$coefficients
mediae
n <- table(course)
n <- 30 
bf1 <- rbind(mediae[3] - mediae[2] - sqrt(2/n * Sp)*qt(1 - alpha / (2*k), DF) ,
             mediae[3] - mediae[2],
             mediae[3]-mediae[2]+ sqrt(2/n * Sp)*qt(1 - alpha / (2*k), DF))

bf2 <- rbind(mediae[3]-mediae[1] - sqrt(2/n * Sp)*qt(1 - alpha / (2*k), DF) ,
             mediae[3]-mediae[1],
             mediae[3]-mediae[1]+ sqrt(2/n * Sp)*qt(1 - alpha / (2*k), DF))
bf3 <- rbind(mediae[2]-mediae[1] - sqrt(2/n * Sp)*qt(1 - alpha / (2*k), DF) ,
             mediae[2]-mediae[1],
             mediae[2]-mediae[1]+ sqrt(2/n * Sp)*qt(1 - alpha / (2*k), DF))
bf1
bf2
bf3
"from group 1 to 3 and 2o to 3 there is significant difference, not any diff
with 2" 
"we know each variance is distributed as:
(ni - 1)si² ~ sigmasq chiquadro(ni-1)"
"Also: (N-g)Sp ~ sigmasq chiquadro(N-G)"
Sp
bf4 <- rbind(c(Sp * DF / qchisq(1 - alpha / (2*k), DF), Sp,
  Sp* DF / qchisq(alpha / (2*k), DF)))
bf4


detach(wait)
#############
### Problem 3
###
sail <- read.table("Sailing.txt")
head(sail)


# define priors
priors <- c(0.8,.2)
attach(sail)
?lda
# priors specified in order of factors
sail.lda <- lda(type ~ as.numeric(water) + as.numeric(sailing.time), prior=priors)
sail.qda <- qda(type ~ as.numeric(water) + as.numeric(sailing.time), prior=priors)

sail2 <- sail[,-3]
i1 <- predict(sail.lda, sail2)
i1 <- i1$class
x11()
plot(sail2, main='Swimmers', xlab='water', ylab='sailing time', pch=20)
points(sail2[which(i1=="seadog"),], col='red', pch=20)
points(sail2[which(i1=="vacationer"),], col='green', pch=20)

legend("topright", legend=levels(as.factor(type)), fill=c('red','green'), cex=.7)

points(sail.lda$means, pch=4,col=c('red','green','blue') , lwd=2, cex=1.5)

x  <- seq(min(sail2[,1]), max(sail2[,1]), length=200)
y  <- seq(min(sail2[,2]), max(sail2[,2]), length=200)
xy <- expand.grid(water=x, sailing.time=y) # need to match

z  <- predict(sail.lda, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - pmax(z[,2])  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
#z2 <- z[,2] - pmax(z[,1])  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    

# Plot the contour line of level (levels=0) of z1, z2, z3: 
# P_i*f_i(x,y)-max{P_j*f_j(x,y)}=0 i.e., boundary between R.i and R.j 
# where j realizes the max.
contour(x, y, matrix(z1,200), levels=0, drawlabels=F, add=T)  
#contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)
#contour(x, y, matrix(z3, 200), levels=0, drawlabels=F, add=T)

x11()
plot(sail2, main='Swimmers', xlab='water', ylab='sailing time', pch=20)
points(sail2[which(i1=="seadog"),], col='red', pch=20)
points(sail2[which(i1=="vacationer"),], col='green', pch=20)
legend("topright", legend=levels(as.factor(type)), fill=c('red','green'), cex=.7)
points(sail.lda$means, pch=4,col=c('red','green','blue') , lwd=2, cex=1.5)
z  <- predict(sail.qda, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - pmax(z[,2])  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
#z2 <- z[,2] - pmax(z[,1])  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    

# Plot the contour line of level (levels=0) of z1, z2, z3: 
# P_i*f_i(x,y)-max{P_j*f_j(x,y)}=0 i.e., boundary between R.i and R.j 
# where j realizes the max.
contour(x, y, matrix(z1,200), levels=0, drawlabels=F, add=T)  



# posteriors
predict(sail.lda, sail[1, -3])$posterior
predict(sail.qda, sail[1, -3])$posterior
lda.pred <- predict(sail.lda)
qda.pred <- predict(sail.qda, sail[1, -3])

## b) we use the APER (we need to adjust with priors)
G <- 2
misc <- table(class.true=type, class.assigned=i1)
APER <- 0
for(g in 1:G)
   APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * priors[g]  

misc <- table(class.true=type, class.assigned=predict(sail.lda, sail2)$class)
APER2 <- 0
for(g in 1:G)
  APER2 <- APER + sum(misc[g,-g])/sum(misc[g,]) * priors[g]  
APER
APER2 # they match, use occam's razor and choose model with less coplexity

z0 <- data.frame(water=35, sailing.time =168)
predict(sail.lda, z0)


###
## Problem 4
##############
lum <- read.table("Lumieres.txt")
head(lum)
attach(lum)
rain <- as.factor(rain)
day2 <- day**2
##
## a)
lm1 <- lm(N.people ~ rain + day + day2 + temperature)
summary(lm1)
x11()
par(c(2,2))
plot(lm1)
"We see a straight line in the qq plot, approximate homoschedasticity, 
no significant leverages that lead to a big cook's distance"
shapiro.test(lm1$residuals)
"Residuals approximately normal, we cannot reject the null hypothesis"

# b)
a <- (c(0,0,1,1,0))
linearHypothesis(lm1, a)
"Reject null hypothesis that it does not have an effect"

# c
lm2 <- lm(N.people ~ rain + day + day2)
summary(lm2)
lm3 <- lm(N.people ~ + day + day2)
summary(lm3)
"Leave as the F test is significant for all coefficients"
## d)
betae <- lm3$coefficients
betae
" test is:"
2*betae[3]*61 + betae[2]
library(car)
a <- c(0, 1, 1*61)
linearHypothesis(lm3, a)
"cannot reject, so we have that they match eatch other"
day3 <- day^2 - 2*61*day
lm4 <- lm(N.people ~ day3 )

z0 <- data.frame(day3 = 58)
predict(lm4, z0, interval = "prediction", level=.95)
