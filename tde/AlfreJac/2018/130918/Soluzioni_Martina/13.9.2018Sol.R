# EXAM 13-09-2018

#################################
######### Exercise 1 ############
#################################

iamg <- read.table('IAMG.txt', header=T)

X1 <- iamg[,1]
X2 <- iamg[,2]
X3 <- iamg[,3]
X <- iamg

# a) Build a confidence region (level 95%) for the mean of X. Characterize the region by reporting 
# its expression, its center, the direction of the axes and the length of the semi-axes.

# Confidence region (centered in x.mean)
# { m \in R^2 s.t. n * (x.mean-m)' %*% (x.cov)^-1 %*% (x.mean-m) < cfr.fisher }
x.mean <- sapply(X,mean)
x.cov <- cov(X)
n <- dim(X)[1]
p <- dim(X)[2]
alpha <- 0.05
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)

# Center:
x.mean

# Directions of the principal axes:
eigen(x.cov/n)$vectors

# Length of the semi-axes of the ellipse:
r <- sqrt(cfr.fisher)
r*sqrt(eigen(x.cov/n)$values) 

# Gaussian assumption on the data
mcshapiro.test(X) # pvalue= 0.22 -> ok

# NOT REQUESTED
# I cannot represents the confidence region in R^3! but I can represents the projection of the 3 directions
# using simultaneous T2 confidenc einterval for example (or using Bonferroni)

T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2

x11()
matplot(1:3,1:3,pch='',ylim=range(X),xlab='Variables',ylab='T2 for a component', 
        main='Simultaneous T2 conf. int. for the components')
for(i in 1:3)segments(i,T2[i,1],i,T2[i,3],lwd=3,col=i)
points(1:3, T2[,2], pch=16, col=1:3)

# with bonferroni correction

k <- p
cfr.t <- qt(1 - alpha/(k*2), n-1)

Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf

# Let's do a plot
x11()
matplot(1:3,1:3,pch='',ylim=range(X),xlab='Variables',ylab='Confidence intervals along a component',main='Confidence intervals')

for(i in 1:3) segments(i,T2[i,1],i,T2[i,3],lwd=2,col='grey35', lty=3)
points(1:3, T2[,1], pch='-', col='grey35')
points(1:3, T2[,3], pch='-', col='grey35')

for(i in 1:3) segments(i,Bf[i,1],i,Bf[i,3],lwd=2,col=i)
points(1:3, Bf[,2], pch=16, col=1:4)
points(1:3, Bf[,1], pch='-', col=1:4)
points(1:3, Bf[,3], pch='-', col=1:4)


# b) Build three T2-simultaneous confidence intervals (level 95%) for: the mean number of 
# registered participants, the mean number of oral presentations and the mean number of no-show.

# I have already done it
T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2

#                 inf  center       sup
# Registered 239.63077 249.44 259.24923
# Talk       115.31917 119.04 122.76083
# No.show     19.46447  23.12  26.77553


# c) Perform a test of level 95% to verify the hypothesis according to which, in mean, only 
# 90% of the registered participants actually show up at IAMG meetings.

# persons who show up: number of registered - number of no show
# H0: 0.9*mu.registered == mu.(reg-no.show) vs H1:H0^C

show.up <- X[,1] - X[,3]

# gaussian assumptions

new.data <- 0.9*X[,1] - show.up 
mean.new <- mean(new.data)
var.new <- var(new.data)

# gaussian assumption
shapiro.test(new.data) # pvalue=0.82 -> ok!

# test for a univariate gaussian population -> t test!
t.test(new.data, alternative = 'two.sided', mu=0, conf.level = 0.95)
# the test show me a pvalue of 0.13 for the alternative hypothesis: true mean is not equal to 0
# So I can accept the null hypothesis and assess that the 90% of the participants show up at IAMG


#################################
######### Exercise 2 ############
#################################

wtimes <- read.table('Waiting.txt', header=T)
attach(wtimes)

# a) Propose a complete ANOVA model for the waiting time as a function of the factors course 
# (starters, main course or dessert) and city (Iasi or Bucarest). Report and verify the assumptions 
# of the model.

# Two-ways ANOVA (complete model, with interaction)
fit.aov <- aov(waiting ~ course + city + course:city)
summary.aov(fit.aov)

# Assumption of the model:
# 1) normality (univariate) in each combination of groups
Ps <- c(shapiro.test(waiting[course=='Starter' & city=='Iasi'])$p,
        shapiro.test(waiting[course=='Starter' & city=='Bucarest'])$p,
        shapiro.test(waiting[course=='Main' & city=='Iasi'])$p,
        shapiro.test(waiting[course=='Main' & city=='Bucarest'])$p,
        shapiro.test(waiting[course=='Dessert' & city=='Iasi'])$p,
        shapiro.test(waiting[course=='Dessert' & city=='Bucarest'])$p)
Ps # ok

# 2) homogeneity of variances
bartlett.test(waiting, course)
bartlett.test(waiting, city)
# ok


# b) Comment on the significance of the factors and of their interaction. If needed, 
# propose a reduced model.

summary.aov(fit.aov)

# In the one at the time tests show in the summary, I can assess that the factor city and the
# interaction between city and course are not significant at any (reasonable) level

# try to remove interaction
fit.aov1 <- aov(waiting ~ course + city)
summary.aov(fit.aov1)

# try to remove city
fit.aov2 <- aov(waiting ~ course)
summary.aov(fit.aov2)

# OK! my reduce model is the one only with course factor
fit.aov.reduce <- aov(waiting ~ course)
summary.aov(fit.aov.reduce)

# c) Build Bonferroni confidence intervals (global level 95%) for the mean differences between 
# the waiting times in the groups identified at point (b), and for the variances of the waiting 
# times within the groups. Comment the results.

N <- dim(wtimes)[1]
g <- length(levels(factor(course)))
DF <- N-g

alpha <- .05
k <- g+1 # +1 because of the interval for the variances

qT <- qt(1-alpha/(2*k), DF)
qCinf <- qchisq(1 - alpha / (2*k), DF)
qCsup <- qchisq(alpha / (2*k), DF)

Spooled <- (t(fit.aov.reduce$res) %*% fit.aov.reduce$res)/DF   
Spooled

m1 <- mean(waiting[which(course=='Starter')])
m2 <- mean(waiting[which(course=='Main')])
m3 <- mean(waiting[which(course=='Dessert')])
medie <- c(m1,m2,m3)

ng <- c(length(which(course=='Starter')),length(which(course=='Main')),length(which(course=='Dessert')))

BF    <- rbind(cbind(inf=medie - sqrt(as.vector(Spooled) / ng) * qT,
                     sup=medie + sqrt(as.vector(Spooled) / ng) * qT),
               c(inf=Spooled * DF / qCinf,
                 sup=Spooled * DF / qCsup))
BF

#       inf      sup
#  77.21149 82.15518  -> mean of Starter
#  79.22815 84.17185  -> mean of Main
#  28.24482 33.18851  -> mean of Dessert
#  44.79582 76.32483  -> variances

# The waiting time for the Dessert is in mean lower then the other.


#################################
######### Exercise 3 ############
#################################

sailing <- read.table('Sailing.txt', header=T)
attach(sailing)

# a) Based on the available features, build two Bayes classifiers, A and B, for the kind of 
# sailboat’ occupation (vacationer, seadog).
# For each classifier report a qualitative plot of the classification regions, and the estimated 
# posterior probability associated with the first observation

# A. the two populations are Gaussian with the same covariance structure -> LDA

group.name <- factor(type, labels=c('seadog','vacationer'))
sail <- sailing[,1:2]

# prior probabilities
p1 <- 0.2
p2 <- 0.8

lda.sail <- lda(sail, group.name, prior = c(p1,p2))
lda.sail

# classification region

x11()
plot(sail, main='sailing type', xlab='water', ylab='sailing.time', pch=20)
points(sail[which(type=='seadog'),], col='red', pch=20)
points(sail[which(type=='vacationer'),], col='green', pch=20)
legend("topright", legend=levels(factor(type)), fill=c('red','green'), cex=.7)

points(lda.sail$means, pch=4,col=c('black','black') , lwd=2, cex=1.5)

x  <- seq(min(sailing[,1]), max(sailing[,1]), length=200)
y  <- seq(min(sailing[,2]), max(sailing[,2]), length=200)
xy <- expand.grid(water=x, sailing.time=y)

z  <- predict(lda.sail, xy)$post  # these are P_i*f_i(x,y)   -> posterior
z1 <- z[,1] - pmax(z[,2])  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- z[,2] - pmax(z[,1])  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}   

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

# estimated of posterior probabilities
Lda.sail <- predict(lda.sail, sail) # useful to have an estimate of the error
names(Lda.sail)
Lda.sail$posterior[1,]
# seadog vacationer 
# 0.95565267 0.04434733 

# B. the two populations are Gaussian with different covariance structures -> QDA

# prior
p1 <- 0.2
p2 <- 0.8

qda.sail <- qda(sail, group.name, prior = c(p1,p2))
qda.sail

# classification region

plot(sail, main='sailing type', xlab='water', ylab='sailing.time', pch=20)
points(sail[which(type=='seadog'),], col='red', pch=20)
points(sail[which(type=='vacationer'),], col='green', pch=20)
legend("topright", legend=levels(factor(type)), fill=c('red','green'), cex=.7)

points(qda.sail$means, pch=4,col=c('black','black') , lwd=2, cex=1.5)

x  <- seq(min(sailing[,1]), max(sailing[,1]), length=200)
y  <- seq(min(sailing[,2]), max(sailing[,2]), length=200)
xy <- expand.grid(water=x, sailing.time=y)

z  <- predict(qda.sail, xy)$post  # these are P_i*f_i(x,y)   -> posterior
z1 <- z[,1] - pmax(z[,2])  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- z[,2] - pmax(z[,1])  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}   

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

# estimated of posterior probabilities
Qda.sail <- predict(qda.sail, sail) # useful to have an estimate of the error
names(Qda.sail)
Qda.sail$posterior[1,]
# seadog vacationer 
# 0.96474706 0.03525294 

# b) Evaluate the performances of the classifiers A and B and identify the best one.
# To evaluate the performances of the two classifiers I compute the APER

# A. 
prior <- c(0.2,0.8)
G <- 2
misc <- table(class.true=group.name, class.assigned=Lda.sail$class)
APER <- 0
for(g in 1:G)
APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior[g]  
APER # 0.063

# B.
prior <- c(0.2,0.8)
G <- 2
misc <- table(class.true=group.name, class.assigned=Qda.sail$class)
APER <- 0
for(g in 1:G)
        APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior[g]  
APER # 0.06


# APER of the LDA is bigger then the QDA one -> QDA is the best one

# c) How would you classify the occupants of a sailboat with daily consumed water 35 l and 
# daily sailing time 168 min?

new.data <- data.frame(water=35 , sailing.time=168)
pred <- predict(qda.sail,new.data)
pred$posterior
# seadog vacationer
# 0.1577169  0.8422831

# I can classify the new observation in the vacationer group


#################################
######### Exercise 4 ############
#################################

lum <- read.table('Lumieres.txt', header=T)
attach(lum)

# Construct the model:

# dummy (=1 if rain, 0 otherwise)
r <- ifelse(rain=='yes',1,0)

fit <- lm(N.people ~ r + day + I(day^2) + temperature)
summary(fit)

# parameters 
params <- data.frame(beta0 = fit$coefficients[1],
                     beta0.rain = fit$coefficients[1] + fit$coefficients[2],
                     beta1 = fit$coefficients[3],
                     beta2 = fit$coefficients[4],
                     beta3 = fit$coefficients[5],
                     sigma=sqrt(sum(residuals(fit)^2)/fit$df))
params

# model assumptions
shapiro.test(fit$residuals) # pvalue=0.8 -> ok

# collinearity:
vif(fit)
# There is collinearity on the two terms days and days^2

# b) Perform a statistical test to verify if the mean number of participants depends 
# significantly on the day of the representation.

# from the summary of the model, looking at one at the time tests (I can do this because the 
# gaussianity assumption on the residuals is verified) I can assess that the parameter
# r is not significant at level 95%

fit.1 <- lm(N.people ~ day + I(day^2) + temperature)
summary(fit.1)

# I can note that the regressor temperature is not significant at this level too
fit.reduce <- lm(N.people ~ day + I(day^2))
summary(fit.reduce)

# c) Based on a statistical test of level 95%, reduce the model and update the parameter estimates.
# the last model

c(fit.reduce$coefficients, sigma=sqrt(sum(residuals(fit.reduce)^2)/fit.reduce$df))

# d) Perform a test to verify if the maximum of the expected number of participants is on the 
# last day of July (d = 61) and, in case, update the estimates of the model parameters.

# derivative: beta1 + 2*beta2*61 = 0

linearHypothesis(fit.reduce, c(0,1,2*61), 0)
# Accept the null hypothesis

# e) Based on the last update of the model parameters, provide a prediction interval (probability 
# 95%) for the number of people participating to the representation on the 28th July 
# (d = 58, temp = 29, no rain).

# beta1 = -2*61*beta2
final.model <- lm (N.people ~ I(-2*61*day + day^2))
summary(final.model)

Z0.new <- data.frame(day=58) 

Pred <- predict(final.model, Z0.new, interval='prediction', level=1-0.05)  
Pred

# fit      lwr      upr
# 309.2318 266.2167 352.2469
