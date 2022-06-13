# EXAM 18-7-2017

################################
######## EXERCISE 1 ############
################################

tourists <- read.table('tourists.txt', header=T)

# a) Perform a Principal Component Analysis of the dataset, by only focusing on the quantitative 
# variables of the dataset; here, evaluate whether it is appropriate to use the original variables 
# or the standardized ones and proceed accordingly. Interpret, when possible, the principal components.
# Report the numerical value of the first 3 loadings and the variance displayed by the data along 
# those directions.

q.tourist <- tourists[,-c(1,2)]

pc.tourists <- princomp(q.tourist, scores=T) #I can compute also the scores
pc.tourists 
summary(pc.tourists)

# standardize variables
tourists.sd <- scale(q.tourist)
tourists.sd <- data.frame(tourists.sd)

pc.tourists.sd <- princomp(tourists.sd, scores=T) 
pc.tourists.sd 
summary(pc.tourists.sd)

# Is better to normalize my data

load.tour <- pc.tourists.sd$loadings
par(mar = c(2,2,2,1), mfrow=c(3,1))
for(i in 1:3)barplot(load.tour[,i], ylim = c(-1, 1), main=paste('Loadings PC ',i,sep=''))

# looking at the loadings, I assess that the first principal component is an overall mean of the variables.
# the second principal component give the contrast between expensive hotel/residence and cheaper one (1-2 
# stars hotel and bed and breakfast). The third principal component take in consideration more the bed and
# breakfast variable.

# Numerical values of first 3 loadings

load.tour[,c(1,2,3)]
#               Comp.1      Comp.2      Comp.3
# Hotel.5stars  0.3409434  0.59123382  0.19309851
# Hotel.4stars  0.3630324  0.30520188 -0.06389270
# Hotel.3stars  0.3682463 -0.07545910 -0.20740553
# Hotel.2stars  0.3519069 -0.42808893 -0.07764867
# Hotel.1star   0.3467628 -0.52252462 -0.17281238
# Residences    0.3580348  0.25373557 -0.30791029
# Bed.Breakfast 0.3379262 -0.17132443  0.86887665
# Rented.flats  0.3604079  0.03967862 -0.17310351

# b) Report (qualitatively) the scatter plot of the data along the first two PCs and describe 
# how to interpret the data clouds in the four quadrants. Use the categorical variables ‘Month’ and 
# ‘Region’ to further interpret the results at point (a).

scores.tourists <- pc.tourists.sd$scores

plot(scores.tourists[,1:2])
abline(h=0, v=0, lty=2, col='grey')



################################
######## EXERCISE 2 ############
################################

horsec <- read.table('horsecolic.txt', header=T)
attach(horsec)

# a) Build 5 Bonferroni confidence intervals (global level 99%) for the mean difference in 
# ‘Rectal.temperature’, ‘Pulse’, ‘Respiratory.rate’ and ‘Packed.cell.volume’ between the horses 
# with pain and without pain. Comment the results and identify the variables along which a significant 
# difference exists. State and verify the appropriate assumptions.

# Prepare the variables

RT.yes <- horsec[which(Pain=='Yes'),1]
RT.no <- horsec[-which(Pain=='Yes'),1]

P.yes <- horsec[which(Pain=='Yes'),2]
P.no <- horsec[-which(Pain=='Yes'),2]

RR.yes <- horsec[which(Pain=='Yes'),3]
RR.no <- horsec[-which(Pain=='Yes'),3]

PC.yes <- horsec[which(Pain=='Yes'),4]
PC.no <- horsec[-which(Pain=='Yes'),4]

t1 <- data.frame(RT.yes,P.yes,RR.yes,PC.yes)
t2 <- data.frame(RT.no,P.no,RR.no,PC.no)

n1 <- dim(t1)[1] # n1=79
n2 <- dim(t2)[1] # n2=221
p  <- dim(t1)[2] # p=4

# assumptions:
mcshapiro.test(t1) # pvalue=0.31 -> ok
mcshapiro.test(t2) # pvalue=0.39 -> ok

# we compute the sample mean, covariance matrices and the matrix Spooled

t1.mean <- sapply(t1,mean)
t2.mean <- sapply(t2,mean)
t1.cov  <-  cov(t1)
t2.cov  <-  cov(t2)
Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2)
# we compare the matrices
list(S1=t1.cov, S2=t2.cov, Spooled=Sp)

# Test H0: mu1 == mu2  vs  H1: mu1 != mu2
# i.e.,
# Test H0: mu1-mu2 == c(0,0,0,0)  vs  H1: mu1-mu2 != c(0,0,0,0)

alpha   <- .01
Spinv   <- solve(Sp)
k <- 4
cfr.t <- qt(1 - alpha/(k*2), n1+n2-1)

Bf <- cbind(inf = t1.mean - t2.mean - cfr.t*sqrt(diag(Sp)*(1/n1+1/n2)),
            center = t1.mean - t2.mean, 
            sup = t1.mean  - t2.mean + cfr.t*sqrt(diag(Sp)*(1/n1+1/n2)))
Bf

#           inf      center        sup
# RT -0.05794228  0.08906925  0.2360808
# P  31.49785036 36.63161178 41.7653732
# RR 16.06946167 19.81537946 23.5612973
# PC -0.46100910  1.83448193  4.1299730
# The variables along which exists a significant difference are PULSE and RESPIRATORY RATE (don't 
# contain zero value)


# b) Based only on the assumptions you deem appropriate, build a classifier for the condition ‘pain’,
# based only on the variables along which the groups display a significant difference according to the
# analysis at point (a). Report the mean within the groups and the prior probabilities estimated from the 
# sample. Report a qualitative plot of the partition induced by the classifier in the space identified by 
# two of the used variables.

pulse.rr <- horsec[,2:3]
yes <- which(Pain=='Yes')
no <- which(Pain=='No')
pain.group <- factor(Pain, labels=c('Yes','No'))

# verify assumptions 1) e 2): 
# 1) normality (univariate) within the groups
mcshapiro.test(pulse.rr[yes,])
mcshapiro.test(pulse.rr[no,])

# 2) equal variance (univariate)
var.test(pulse.rr[yes,1],pulse.rr[no,1]) 
var.test(pulse.rr[yes,2],pulse.rr[no,2])
# I can't accept this assumption! I use QDA

nA <- length(yes)
nB <- length(no)
n  <- nA + nB

# Prior probabilities (estimated from the data, no prior knowledge)
qda.horse <- qda(pulse.rr, pain.group)
qda.horse

plot(pulse.rr, main='pulse and rr', xlab='p', ylab='rr', pch=20)
points(pulse.rr[yes,], col='red', pch=20)
points(pulse.rr[no,], col='green', pch=20)
points(qda.horse$means, col=c('black','black'), pch=4, lwd=2, cex=1.5)

x  <- seq(min(pulse.rr[,1]), max(pulse.rr[,1]), length=200)
y  <- seq(min(pulse.rr[,2]), max(pulse.rr[,2]), length=200)
xy <- expand.grid(p=x, rr=y)

z  <- predict(qda.horse, xy)$post  
z1 <- z[,1] - z[,2]    

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)

# c) Estimate the APER of classifier.
Qda.horse <- predict(qda.horse, pulse.rr)
errorsCV <- (Qda.horse$class != pain.group) # I can use this method because I use for prior probabilities
                                            # the empirical frequencies
APER   <- sum(errorsCV)/length(pain.group)
APER # 0.053


################################
######## EXERCISE 3 ############
################################

castle <- read.table('castle.txt', header=T)

# a) Perform a statistical test to verify if the centre of the distribution is located in the 
# centre of Aosta (Lat = 45.733, Long = 7.333). Verify the assumptions of the test.

# H0: mu == c(45.733, 7.333) vs H1: mu != c(45.733, 7.333)

# Gaussian hypothesis
mcshapiro.test(castle) #pvalue= 0.4 -> OK

n <- dim(castle)[1]
p <- dim(castle)[2]

x.mean <- sapply(castle, mean)
x.cov <- cov(castle)
x.invcov <- solve(x.cov)

alpha <- 0.05
mu0 <- c(45.733,7.333)

# T2 Statistics
x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
# Radius of the ellipsoid
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)  # qf-> compute the quantile of F distribution 
# Test: 
x.T2 < cfr.fisher   # TRUE -> no statistical evidence to reject H0 at level alpha

# Compute the p-value 
P <- 1-pf(x.T2*(n-p)/((n-1)*p), p, n-p)
P # = 0.28

# So I accept the null hypothesis

# b) Consistently with the results at point (a), estimate an elliptical region that contains 95% 
# of the castles. Report: the analytical expression of the region, its centre, the direction and 
# the length of the principal axes of the ellipse. Report a qualitative plot of the region.

plot(castle, asp = 1)

# Center:
x.mean 
# Lat      Long 
# 45.731404  7.334415 

# Radius of the ellipse:
r <- sqrt(qchisq(0.95,2)) # It's a PREDICTION region! 
r # 2.447747

# Length of the semi-axes:
r*sqrt(eigen(x.cov)$values)
# 0.09661096 0.01471144

ellipse(x.mean, x.cov, r, col = 'red', lty = 2, lwd=2)

################################
######## EXERCISE 4 ############
################################

albatross <- read.table('albatross.txt', header=T)
attach(albatross)

# Dummy variable:
upwind <- ifelse(wind=='upwind',1,0)
va.2 <- I(Va^2)
vi.2 <- I(Vi^2)

# Model:
model <- lm(distance ~ va.2 + vi.2 + upwind + va.2:upwind + vi.2:upwind)
summary(model)

# a) Estimate the 7 parameters of the model (report αg , βg , γg for g = 1, 2 and σ) and verify
# its assumptions.

# coefficients (d=downwind, u=upwind)
coeffs <- data.frame(alpha.d = model$coefficients[1],
                     alpha.u = model$coefficients[1] + model$coefficients[4],
                     beta.d = model$coefficients[2],
                     beta.u = model$coefficients[2] + model$coefficients[5],
                     gamma.d = model$coefficients[3],
                     gamma.u = model$coefficients[3] + model$coefficients[6],
                     sigma=sqrt(sum(residuals(model)^2)/model$df))
# alpha.d  alpha.u     beta.d     beta.u    gamma.d      gamma.u    sigma
# 3.288063 4.651341 0.01970956 0.01114171 -0.0175382 -0.009772231 1.799546

# gaussian hypothesis on the residuals
shapiro.test(model$residuals) #pvalue = 0.5952 -> OK
# omoschedasticity
plot(model$residuals) # no patterns! ok

# b) Based on appropriate test(s), reduce the model.
# I have verified the gaussian assumption on the residuals of the model, so I can observe the
# one at the time tests in the summary of the model:
# the regressor "upwind" is not significant! I can remove it:

model <- lm(distance ~ va.2 + vi.2 + va.2:upwind + vi.2:upwind)
summary(model)
# all the others are significant at level 10%

# c) Using model (b), test the hypothesis according to which γg = −βg , for g = 1, 2. Possibly
# propose a constrained model and estimate its parameters.
# H0: gamma.g + beta.g == 0

linearHypothesis(model, rbind(c(0,0,0,1,1),
                              c(0,1,1,0,0)),c(0,0))
# I have statistical evidence to accept null hypothesis!

model.reduce <- lm(distance ~ I(va.2 - vi.2) + I(va.2 - vi.2):upwind)
summary(model.reduce)


# Wilbur arbatross is giving Bianca and Bernie a lift to Australia. Do you deem the landing 
# of Bianca and Bernie to be safe in case of upwind or downwind wind, on a 17 m long runway, 
# if Wilbur is approaching with Va = 35 km/h and Vi = 25 km/h? Base your answer on model (c) and 
# on two intervals of global level 99% for Wilbur’s landing distance.

alpha <- 0.01/2 # I use bonferroni correction
z_upwind <- data.frame(va.2 = 35^2, vi.2 = 25^2, upwind = 0)
z_downwind <- data.frame(va.2 = 35^2, vi.2 = 25^2, upwind = 1)
predict(model.reduce, z_upwind, interval = 'prediction',level = 1-alpha)
predict(model.reduce, z_downwind, interval = 'prediction',level = 1-alpha)

# upwind:
#     fit     lwr      upr
# 16.25067 10.9296 21.57175

# downwind:
#     fit      lwr      upr
# 11.59718 6.277036 16.91733

# Is safer with downwind! 





