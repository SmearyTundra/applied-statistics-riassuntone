#######
# REGRESSION
######

library(MASS)
library(car)
library(rgl)


library(MASS)
library(car)
library(rgl)
library(glmnet)
library(ISLR)
library(leaps)
library(tree)


##### Basic regression + linear hyp + pred ####

airfoil <- read.table('airfoil.txt', header=T)
head(airfoil)
airfoil$velocity <- factor(airfoil$velocity)

### Model:
### distance = beta_0 + beta_1 * speed + beta_2 * speed^2 + Eps
### (linear in the parameters!)

### Assumptions:
## 1) Parameter estimation: E(Eps) = 0  and  Var(Eps) = sigma^2 
## 2) Inference:            Eps ~ N(0, sigma^2)

### Assumptions Parameter estimation: E(Eps) = 0  and  Var(Eps) = sigma^2 

#regression
fit1 <- lm(sound ~ frequency*velocity, data = airfoil)   #I(velocity^2)
summary(fit1)
summary(fit1)$sigma^2 #sum(residuals(fit1)^2)/fit1$df  # estimate of sigma^2
fit1$coefficients
shapiro.test(fit1$residuals)
vif(fit1) # high => highly collinear with the other variables in the model
par(mfrow=c(2,2))
plot(fit1)

##### Assumption Inference: Eps ~ N(0, sigma^2)

#linear hyp
linearHypothesis(fit1, rbind(c(0,1,0,0), c(0,0,0,1)),c(0,0))
linearHypothesis(fit1, rbind(c(0,0,1,0), c(0,0,0,1)),c(0,0))
linearHypothesis(fit1, rbind(c(0,0,0,1)),c(0))

fit2 <- lm(sound ~ frequency + velocity, data = airfoil)
summary(fit2)
shapiro.test(fit2$residuals)
vif(fit2)
par(mfrow=c(2,2))
plot(fit2)

#prediction & confidence interval 
new <- data.frame(frequency = 15000, velocity = 'H' )
# or new <- data.frame(name = seq(min(original_data),max(original_data),len=100))
answer <- predict(fit2, newdata = new, interval = 'confidence') #level = 1-alpha
                                                  #'prediction'
                                                  #'
#confidence interval for mean and variance
Z0   <- data.frame(profondita = 200, pozzo = '3')
ICBmean <- predict(fit2, Z0, interval='confidence',level=1-alpha/k) 

e <- residuals(fit2)
ICBvar <- data.frame(L=t(e)%*%e/qchisq(1-alpha/(2*k),n-(r+1)),
                     U=t(e)%*%e/qchisq(alpha/(2*k),n-(r+1)))
ICBvar

#max of the values & prediction of a data already in the dataset
t.max <- which.max(result4$fitted.values)
max <- predict(result4, df, interval = 'confidence', level = 0.99)[t.max,]

# Bonferroni confidence intervals on beta
p <- length(fm$coefficients)-1  # number of tested coefficients
confint(fm, level= 1-0.05/p)[2:3,]  # Bonferroni correction!
# Note: `confint()` returns the confidence intervals one-at-a-time; to have a global level 95% we need to include a correction







# COLLINEARITY


##### PCA regression
speed.pc <- princomp(cbind(speed1,speed2), scores=TRUE)
summary(speed.pc)
speed.pc$loadings

sp1.pc <- speed.pc$scores[,1]
sp2.pc <- speed.pc$scores[,2]

fm.pc <- lm(distance ~ sp1.pc + sp2.pc) # fit the PCs
summary(fm.pc) 

##### Ridge regression ####

lambda <- .5
fit.ridge <- lm.ridge(distance ~ speed1 + speed2, lambda = lambda)
names(fit.ridge) # properties we can access
yhat.lm <- cbind(rep(1,n), speed1, speed2) %*% coef(fm)  # LM fitted values
yhat.r  <- cbind(rep(1,n), speed1, speed2) %*% coef(fit.ridge) # ridge fitted values

# test many lambdas
lambda.c <- seq(0,10,0.01)
fit.ridge <- lm.ridge(distance ~ speed1 + speed2, lambda = lambda.c) # lambda is a sequence here!
select(fit.ridge) # choose the best one


##### Lasso regression ####

# Build the matrix of predictors
x <- model.matrix(tox~., data = toxicity)[,-1]
# Build the vector of response
y <- toxicity$tox

# Let's set a grid of candidate lambda's for the estimate
lambda.grid <- seq(0.01,1,length=100)
fit.lasso <- glmnet(x,y, lambda = lambda.grid, alpha = 1)
# alpha = 1 -> Lasso (default)
# alpha = 0 -> Ridge

# Plot coefficients value versus lambda
plot(fit.lasso,xvar='lambda',label=TRUE, col = rainbow(dim(x)[2]))
legend('topright', dimnames(x)[[2]], col =  rainbow(dim(x)[2]), lty=1, cex=1)

# Let's set lambda via cross validation
cv.lasso <- cv.glmnet(x,y,lambda=lambda.grid, nfolds = 10) # default: 10-fold CV

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso
plot(cv.lasso)
abline(v=log(bestlam.lasso), lty=1)
# solid line: optimal lambda
# dashed line: biggest lambda not to far from the optimal one (1 std from the optimal)

# Get the coefficients for the optimal lambda
coef.lasso <- predict(fit.lasso, s=bestlam.lasso, type = 'coefficients')
coef.lasso


# another way to do model selection
library(leaps)
help(regsubsets)

regfit.full <- regsubsets(Salary ~ ., data=Hitters,
                          nvmax=19, # max number of subsets
                          method="forward-backward-exhaustive")
summary(regfit.full)
reg.summary$which
reg.summary$rsq   # r-squared
reg.summary$adjr2 # adjusted r-squared
reg.summary$rss   # residual sum of squares

best_model_id <- which.max(reg.summary$adjr2)
coef(regfit.full, best_model_id)
