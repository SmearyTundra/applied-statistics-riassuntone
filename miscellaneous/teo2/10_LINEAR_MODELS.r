### TOPICS:
### Linear models

library(MASS)
library(car)
library(rgl)


############### LINEAR MODELS #############
######------------------------------------
fm <- lm(distance ~ speed1 + speed2)

summary(fm) 
# Note: collinearity

fitted(fm)        # y hat
residuals(fm)     # eps hat

coefficients(fm)  # beta_i
vcov(fm)          # cov(beta_i)

fm$rank # order of the model [r+1]
fm$df   # degrees of freedom of the residuals [n-(r+1)]

hatvalues(fm) # h_ii
rstandard(fm) # standardized residuals

summary(fm)$sigma^2  # estimate of sigma^2



##### Diagnostic
##### Assumption: Eps ~ N(0, sigma^2)
#####-------------------------------------------
par(mfrow=c(2,2))
plot(fm)

vif(fm)

shapiro.test(residuals(fm))




##### Linear Hypotheses Testing
#####------------------------------------------
### Test (Fisher):
# H0: (beta1, beta2) == (0, 0) vs H1: (beta1, beta2) != (0, 0)
linearHypothesis(fm, rbind(c(0,1,0), 
						   c(0,0,1)), c(0,0)) 




##### Confidence region for betas
#####-----------------------------------------
p <- 2  # number of tested coefficients
r <- 2  # number of regressors

# Confidence region:
# center: point estimate
c(coefficients(fm)[2], coefficients(fm)[3])
# Direction of the axes?
eigen(vcov(fm)[2:3,2:3])$vectors

plot(coefficients(fm)[2], coefficients(fm)[3], xlim = c(-6,6), ylim = c(-6,6), asp=1, xlab='beta1', ylab='beta2')
ellipse(coefficients(fm)[2:3], vcov(fm)[2:3,2:3], sqrt(p*qf(1-0.05,p,n-(r+1))))
abline(v=0)
abline(h=0)
# Note: collinearity!




##### Bonferroni intervals on betas
#####-----------------------------------------
# Bonferroni intervals (level 95%)

# here testing only coef 2 and 3, not 1 (intercept)

Bf <- rbind(
  beta1=c(coefficients(fm)[2]-sqrt(vcov(fm)[2,2])*qt(1-0.05/(2*p), n-(r+1)),
          coefficients(fm)[2],
          coefficients(fm)[2]+sqrt(vcov(fm)[2,2])*qt(1-0.05/(2*p), n-(r+1))),
  beta2=c(coefficients(fm)[3]-sqrt(vcov(fm)[3,3])*qt(1-0.05/(2*p), n-(r+1)),
  		  coefficients(fm)[3],
          coefficients(fm)[3]+sqrt(vcov(fm)[3,3])*qt(1-0.05/(2*p), n-(r+1)))
)
Bf

# ALTERNATIVE (only for intervals on beta)
p <- length(fm$coefficients)-1  # number of tested coefficients
confint(fm, level= 1-0.05/p)[2:3,]  # Bonferroni correction!

# Note: confint() returns the confidence intervals one-at-a-time;
# to have a global level 95% we need to include a correction



##### Bonferroni intervals on sigma^2
#####-----------------------------------------
# Bonferroni intervals (level 1-alpha)
e <- residuals(fm)
ICBvar <- data.frame(L=t(e)%*%e/qchisq(1-alpha/(2*k),n-(r+1)),
                     U=t(e)%*%e/qchisq(alpha/(2*k),n-(r+1)))
ICBvar



##### Confidence intervals for the mean
##### & prediction (new obs)
##### Assumption: Eps ~ N(0, sigma^2)
#####---------------------------------
# Command predict()

Z0.new <- data.frame(speed1=10, speed2=10^2)

# Conf. int. for the mean
Conf <- predict(fm, Z0.new, interval='confidence', level=1-0.05)  
Conf
# Pred. int. for a new obs
Pred <- predict(fm, Z0.new, interval='prediction', level=1-0.05)  
Pred


# MAX OF THE VALUES & PREDICTION OF A DATA ALREADY IN THE DATASET
t.max <- which.max(fm$fitted.values)
max <- predict(fm, data, interval = 'confidence', level = 0.99)[t.max,]

