piadeina <- read.table("piadeina.txt", header = T)
#create dummy variable:
D <- as.factor(piadeina$Day.of.Week)
attach(piadeina)
g <- lm(Sales ~ D + Bread.Sold + Wraps.Sold + Sandwich.Sold + Focaccia.Sold + Piadina.Sold + Chips.Sold + Juices.Sold + Total.Soda.and.Coffee.Sold + Max.Daily.Temperature)
summary(g)
betas <- g$coefficients
betas[2:5] <- betas[2:5] + betas[1]
names(betas)[1] <- "DFri"
sigma <- summary(g)$sigma

#verification of the assumptions:
shapiro.test(g$residuals)
x11()
qqnorm(g$residuals)
qqline(g$residuals)
x11()
plot(g$fitted.values, g$residuals/sigma)
x11()
boxplot(g$residuals/sigma ~ D)
#all good



#we now fit Lasso regression with penalization coefficient 5
library(glmnet)
x <- model.matrix(Sales ~ D + Bread.Sold + Wraps.Sold + Sandwich.Sold + Focaccia.Sold + Piadina.Sold + Chips.Sold + Juices.Sold + Total.Soda.and.Coffee.Sold + Max.Daily.Temperature)[,-1]
y <- Sales

fit.lasso <- glmnet(x,y, lambda = 5) # default: alpha=1 -> lasso 
#we get the coefficients:
coef.lasso <- predict(fit.lasso, type = 'coefficients')
coef.lasso 

#now for cross validation
lambda.grid <- seq(length = 1000,0,100)

cv.lasso <- cv.glmnet(x,y,lambda=lambda.grid) # default: 10-fold CV

plot(cv.lasso)
bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso
#in this case 6.606607, but it depends on the samples selected by 10-fold
fit.lasso.opt <- glmnet(x,y, lambda = bestlam.lasso) # default: alpha=1 -> lasso 
#we get the coefficients:
coef.lasso.opt <- predict(fit.lasso.opt, type = 'coefficients')
coef.lasso.opt
