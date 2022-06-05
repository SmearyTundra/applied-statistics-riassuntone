#######
# REGRESSION
######

##### Basic regression + linear hyp + pred ####

airfoil <- read.table('airfoil.txt', header=T)
head(airfoil)
airfoil$velocity <- factor(airfoil$velocity)

### Assumptions Parameter estimation: E(Eps) = 0  and  Var(Eps) = sigma^2 

#regression
fit1 <- lm(sound ~ frequency*velocity, data = airfoil)   #I(velocity^2)
summary(fit1)
shapiro.test(fit1$residuals)
vif(fit1)
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

##### Lasso regression ####

# Build the matrix of predictors
x <- model.matrix(tox~., data = toxicity)[,-1]
# Build the vector of response
y <- toxicity$tox

# Let's set a grid of candidate lambda's for the estimate
lambda.grid <- seq(0.01,1,length=100)
fit.lasso <- glmnet(x,y, lambda = lambda.grid) # default: alpha=1 -> lasso 
                                               #if alpha=0 -> ridge regression

plot(fit.lasso,xvar='lambda',label=TRUE, col = rainbow(dim(x)[2]))
legend('topright', dimnames(x)[[2]], col =  rainbow(dim(x)[2]), lty=1, cex=1)
# Let's set lambda via cross validation
cv.lasso <- cv.glmnet(x,y,lambda=lambda.grid) # default: 10-fold CV

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso
plot(cv.lasso)
abline(v=log(bestlam.lasso), lty=1)

# Get the coefficients for the optimal lambda
coef.lasso <- predict(fit.lasso, s=bestlam.lasso, type = 'coefficients')
coef.lasso 