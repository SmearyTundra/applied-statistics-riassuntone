# EXAM 13/09/2018 ex 4

setwd("D:/Dati/Ele/poli/anno4/applied statistics/OldExams/tde")
load("D:/Dati/Ele/poli/anno4/applied statistics/R lab/LAB_5/mcshapiro.test.RData")

library(MASS)
library(carData)
library(car)
library(rgl)

# The file Lumieres.txt reports the data on the participation to Rendez-vous,
# the show of sounds and lights which takes place every evening between the
# 1st June and the 31st August in the main square of Nancy (France). The
# dataset refers to the years 2016 to 2018, and reports:
# - the number of participants (n),
# - the day of the representation (d, d=1 on the 1st June,...,d=92 on the 31st 
#   August),
# - the temperature recorded that evening at 10 pm (temp) and 
# - the weather conditions (rain: yes/no). 
# Consider the following model
# n.ig = beta0.g + beta1*d.i + beta2*d.i^2 + beta3*temp + eps
# where g???{1,2} indicates the group according to the weather conditions
# (g = 1 for no rain, g = 2 for rain) and eps ??? N(0, sigma^2).

lum <- read.table('Lumieres.txt', header=T)
head(lum)
attach(lum)

# a) Estimate the 6 parameters of the model. Verify the model assumptions.

# We have a categorical variable -> g=2 groups -> g-1=1 dummy variable
dummy <- ifelse(rain=='yes', 1, 0) # Dummy for rain

# Model (with dummy variable):
# n.i = b0 + b1*dummy + b2*d.i + b3*d.i^2 + b4*temp + eps
# where
# beta0.1=b0
# beta0.2=b0+b1
# beta1=b2
# beta2=b3
# beta3=b4

fit <- lm(N.people ~ dummy + day + I(day^2) + temperature, data=lum)
summary(fit)
coeffs <- data.frame(beta0.1=coef(fit)[1],
                     beta0.2=coef(fit)[1]+coef(fit)[2],
                     beta1=coef(fit)[3],
                     beta2=coef(fit)[4],
                     beta3=coef(fit)[5],
                     sigma= sqrt(sum(residuals(fit)^2)/fit$df.residual))
coeffs
# beta0.1   beta0.2   beta1     beta2        beta3       sigma
# 174.1094  168.8532  5.007955  -0.04295472  -0.3187537  21.1601

# Assumptions:
# Gaussianity of the residuals:
shapiro.test(residuals(fit)) # pvalue=0.8108 -> OK
# Homoschedasticity
par(mfrow=c(2,2))
plot(fit) # No patterns in the residuals plot -> OK
# Independence
vif(fit) # There is collinearity between the 2 terms days and days^2

# __________________________________________________________________________

# b) Perform a statistical test to verify if the mean number of participants
#    depends significantly on the day of the representation.

# Test: H0: (b2,b3)==(0,0)   vs   H1: (b2,b3)!=(0,0)
linearHypothesis(fit,
                 rbind(c(0,0,1,0,0),
                       c(0,0,0,1,0)),
                 c(0,0))
# Pvalue<2.2e-16 -> Very small -> Reject H0 for any alpha -> The mean number of 
# participants depends significantly on the day of the representation.
# We can see this also from the summary(fit)

# _________________________________________________________________________

# c) Based on a statistical test of level 95%, reduce the model and update the
#    parameter estimates.

# Looking at the pvalues of the one-at-a-time tests in the summary(fit), we see
# that some of the variables are not significant (they have big pvalues), so 
# we can try to reduce the model.

# First we remove temperature (biggest pvalue):
fit2 <- lm(N.people ~ dummy + day + I(day^2), data=lum)
summary(fit2)
# Variable dummy is not significant at 95% (pvalue=0.433)

# Remove dummy
fit3 <- lm(N.people ~ day + I(day^2), data=lum)
summary(fit3)
# All pvalues are small -> All variables are significant

# So we obtained the reduced model:
# n.i = b0 + b1*d.i + b2*d.i^2 + eps
# where
# beta0.1=b0
# beta0.2=b0
# beta1=b1
# beta2=b2

coeffs.red <- data.frame(beta0.1=coef(fit3)[1],
                         beta0.2=coef(fit3)[1],
                         beta1=coef(fit3)[2],
                         beta2=coef(fit3)[3],
                         sigma= sqrt(sum(residuals(fit3)^2)/fit3$df.residual))
coeffs.red
# beta0.1   beta0.2    beta1      beta2         sigma
# 164.2542  164.2542   5.024694   -0.04296185   20.90942

# ______________________________________________________________________

# d) Perform a test to verify if the maximum of the expected number of
#    participants is on the last day of July (d = 61) and, in case, update the
#    estimates of the model parameters.

# We have to make a computation with the maximum -> We need to find the
# derivative of our model

# Reduced model: n.i = b0 + b1*d.i + b2*d.i^2 + eps
# Derivative: b1 + 2*b2*d.i

# We impose the derivative equal to zero to obtain the equation for the maximum:
# b1 + 2**b2*d.i = 0

# Test: H0: b1 + 2*b2*61 == 0   vs   H1: b1 + 2*b2*61 != 0
linearHypothesis(fit3,
                 c(0,1,2*61),
                 0)
# Pvalue=0.2073 -> Quite big -> Do not reject H0 -> The relation is correct ->
# The maximum of the expected number of participants is on the last day of July 

# Thus we can propose a new constrained model, substituting this relationship
# in the previous expression:

# Model (with constraint):
# n.i = b0 - 2*61*b2*d.i + b2*d.i^2 + eps
#     = b0 + (d.i^2-2*61*d.i)*b2 + eps
# with eps ~ N(0, sigma^2)
fit4 <- lm(N.people ~ I(day^2-2*61*day), data=lum)
summary(fit4)
coeffs.con <- data.frame(b0=coef(fit4)[1],
                         b2=coef(fit4)[2],
                         sigma= sqrt(sum(residuals(fit4)^2)/fit4$df.residual))
coeffs.con
#       b0          b2    sigma
# 166.1655 -0.03854158 21.03339

# ____________________________________________________________________________

# e) Based on the last update of the model parameters, provide a prediction
#    interval (probability 95%) for the number of people participating to the
#    representation on the 28th July (d = 58, temp = 29, no rain).

alpha <- 0.05
z0.new <- data.frame(day=58)
IP <- predict(fit4, newdata = z0.new, interval = 'prediction', level = 1-alpha)
IP
# fit      lwr      upr
# 309.2318 266.2167 352.2469
 