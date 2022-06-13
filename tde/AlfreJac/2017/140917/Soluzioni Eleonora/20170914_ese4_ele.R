# EXAM 14/9/2017 ex 3

setwd("D:/Dati/Ele/poli/anno4/applied statistics/OldExams/tde")
load("D:/Dati/Ele/poli/anno4/applied statistics/R lab/LAB_5/mcshapiro.test.RData")

library(MASS)
library(car)
library(rgl)

# The Venice lagoon is affected by important tide phenomena, and particularly
# by the so-called acqua alta. Consider the following model for the
# sea level H [cm] at Punta della Salute (in Venice) at 17:00 of each day:
# H = Beta.0 + Beta.1*sin(2*pi/28*t) + Beta.2*sin(pi/365*(t-t0)) + Beta.3*t + eps
# where:
# t belonging to [1, 365] is the day of the year,
# t0 = 82 indicates the vernal equinox,
# eps ~  N(0, sigma^2).
# Interpret:
# - the first term as due to the effect of the moon motion (astronomical tide),
# - the second term as associated with seasonal effects,
# - the third as due to the global increase of the sea level.
# The file tide.txt reports the sea level measured at 17:00 on 203 days of 2015.

tide <- read.table('tide.txt', header=T)
head(tide)
attach(tide)

# a) Estimate the five parameters of the model. Report the estimates of Beta.i, 
#    i in {0,1,2,3}, and of sigma.

# Model:
# H = Beta.0 + Beta.1*sin(2*pi/28*t) + Beta.2*sin(pi/365*(t-t0)) + Beta.3*t + eps
fit <- lm(h ~ I(sin(2*pi/28*t)) + I(2*sin(pi/365*(t-82))) + t, data=tide)
summary(fit)

# Assumptions/Diagnostics
shapiro.test(residuals(fit)) # pvalue=0.2402 -> Gaussianity of residuals: OK
x11()
par(mfrow=c(2,2))
plot(fit) # All ok
# Gaussianity and homoschedasticity of residuals: OK

coefficients(fit)
sigma <- sqrt(sum(residuals(fit)^2)/fit$df.residual)
sigma
# Beta.0 = (Intercept) = 66.49527511
# Beta.1 = I(sin(2 * pi/28 * t)) = 19.21718913
# Beta.2 = I(2 * sin(pi/365 * (t - 82))) = 0.84866532 
# Beta.3 = t = 0.02371382 
# sigma = 8.904655

# _______________________________________________________________________

# b) Having introduced and verified the appropriate assumptions, perform two
#    statistical tests to verify if
#    - the mean sea level is influenced by the periodic components;
#    - the mean sea level depends on the global increase of the sea level.

# Test 1:
# H0: (Beta.1, Beta.2)==(0, 0)   vs   H1: (Beta.1, Beta.2)!=(0, 0)
linearHypothesis(fit,
                 rbind(c(0,1,0,0),
                       c(0,0,1,0)),
                 c(0,0))
# Pvalue<2.2e-16 -> Reject H0 for any alpha -> The periodic terms are relevant

# However, looking at the one-at-a-time pvalues of summary(fit) (we apply 
# Bonferroni correction), we see that only one of them is significant:
# I(sin(2 * pi/28 * t)): pvalue<2e-16 -> Significant
# I(2 * sin(pi/365 * (t - 82))): pvalue=0.5130/2=0.2565 -> Not significant
# -> We can remove the not significant term.

# Test 2:
# H0: Beta.0==0   vs   H1: Beta.0!=0
linearHypothesis(fit,
                 c(1,0,0,0),
                 0)
# Pvalue<2.2e-16 -> Reject H0 for any alpha -> The global increase of the sea
# level is significant

# ____________________________________________________________________________

# c) Based on point (b), propose a reduced model and estimate its parameters.

# We can reduce the model, removing the not significant term
# Reduced model:
# H = A.0 + A.1*sin(2*pi/28*t) + A.2*t + eps
fit2 <- lm(h ~ I(sin(2*pi/28*t)) + t, data=tide)
summary(fit2)

coefficients(fit2)
sigma2 <- sqrt(sum(residuals(fit2)^2)/fit2$df.residual)
sigma2
# A.0 = (Intercept) = 66.01383194
# A.1 = I(sin(2 * pi/28 * t)) = 19.18870836
# A.2 = t = 0.03082913 
# sigma2 = 8.891944

# _________________________________________________________________________

# d) Based on model (c), provide two prediction intervals (global level 90%)
#    for the sea level at 17:00 of 20th September 2017 (day 263 of 2017) and of
#    1st December 2017 (day 335 of 2017). Comment the results knowing that in
#    Venice high-water is expected whenever the sea level is higher than 90 cm.

z0.new <- data.frame(t=c(263, 335))
alpha <- 0.10
k <- 2 # Correction to make the two predictions hold simultaneously

IP <- predict(fit2, newdata = z0.new, interval = 'prediction', level = 1-alpha/k)
IP
# fit      lwr      upr
# 86.08586 68.45025 103.7215   20th September 2017
# 72.07170 54.40830  89.7351   1st December 2017

# With overall confidence 90%, we won't have high water on 1st December 2017,
# while we may have it on 20th September 2017

# If instead we do not want the two predictions to hold simultaneously:
# 20th September 2017
z0.new.1 <- data.frame(t=263)
IP.1 <- predict(fit2, newdata = z0.new.1, interval = 'prediction', level = 1-alpha)
IP.1
#      fit      lwr     upr
# 86.08586 71.30669 100.865

# 1st December 2017
z0.new.2 <- data.frame(t=335)
IP.2 <- predict(fit2, newdata = z0.new.2, interval = 'prediction', level = 1-alpha)
IP.2
# fit      lwr      upr
# 72.0717 57.26924 86.87416

# Conclusions do not change