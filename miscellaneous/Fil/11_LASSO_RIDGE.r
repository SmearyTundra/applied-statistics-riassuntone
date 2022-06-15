### TOPICS:
### Linear models - LASSO - RIDGE

library(MASS)
library(car)
library(rgl)

library(glmnet)


# Variance inflation factor
vif(fm)

### A possible solution to collinearity: PCA
speed.pc <- princomp(cbind(speed1,speed2), scores=TRUE)
summary(speed.pc)
speed.pc$load

sp1.pc <- speed.pc$scores[,1]
sp2.pc <- speed.pc$scores[,2]

# Now we estimate the model by inserting the PCs instead of the 
# original regressors 
# Model: y = b0 + b1*PC1+ b2*PC2 + eps, eps~N(0,sigma^2)
fm.pc <- lm(distance ~ sp1.pc + sp2.pc)

summary(fm.pc) 




### Another possible solution: RIDGE
# Fix lambda
lambda <- .5
fit.ridge <- lm.ridge(distance ~ speed1 + speed2, lambda = lambda)
# Note: to fit the model, R automatically centers X and Y 
# with respect to their mean.
coef.ridge <- coef(fit.ridge)
yhat.lm <- cbind(rep(1,n), speed1, speed2)%*%coef(fm)  # LM fitted values
yhat.r <- cbind(rep(1,n), speed1, speed2)%*%coef.ridge # ridge fitted values


# Choice of the optimal lambda, e.g., via cross-validation
select(fit.ridge)

# or
lambda.opt <- lambda.c[which.min(fit.ridge$GCV)]
lambda.opt





### Another possible solution: LASSO
# Build the matrix of predictors
x <- model.matrix(distance~speed1+speed2)[,-1]
# Build the vector of response
y <- distance

# Let's set a grid of candidate lambda's for the estimate
lambda.grid <- 10^seq(5,-3,length=100)
fit.lasso <- glmnet(x,y, lambda = lambda.grid) # default: alpha=1 -> lasso 
# [note: if alpha=0 -> ridge regression]

plot(fit.lasso,xvar='lambda',label=TRUE, col = rainbow(dim(x)[2]))
legend('topright', dimnames(x)[[2]], col =  rainbow(dim(x)[2]), lty=1, cex=1)

# Let's set lambda via cross validation
cv.lasso <- cv.glmnet(x,y,lambda=lambda.grid) # default: 10-fold CV

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

plot(cv.lasso)
abline(v=log(bestlam.lasso), lty=1)

# Get the coefficients for the optimal lambda
coef.lasso <- predict(fit.lasso, s=bestlam.lasso, type = 'coefficients')[1:3,]
coef.lasso 