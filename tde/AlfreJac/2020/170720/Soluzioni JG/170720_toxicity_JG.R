toxicity <- read.table("toxicity.txt", header = T)
attach(toxicity)
#model:
#toxicity = beta_0 + sum{i = 1 : 6} (beta_i*Ci)
g <- lm(tox ~ C1 + C2 + C3 + C4 + C5 + C6)
summary(g)

#verification of the assumptions:
shapiro.test(g$residuals)
x11()
qqnorm(g$residuals)
qqline(g$residuals)
#we assume gaussianity
x11()
plot(g$fitted.values, g$residuals/summary(g)$sigma) #few bad residuals, but still acceptable

x11()
layout(rbind(c(1,2,3),c(4,5,6)))
plot(C1, g$residuals/summary(g)$sigma)
plot(C2, g$residuals/summary(g)$sigma)
plot(C3, g$residuals/summary(g)$sigma)
plot(C4, g$residuals/summary(g)$sigma)
plot(C5, g$residuals/summary(g)$sigma)
plot(C6, g$residuals/summary(g)$sigma)
#no patterns

betas<- g$coefficients
s <- summary(g)$sigma

z_new <- data.frame(C1=100, C2=0.7, C3=2, C4=4, C5=1.4, C6=3)
alpha <- 0.05
#we are interested in a PREDICTION interval
predict(g, z_new, interval = "prediction", level = 1 - alpha)



library(glmnet)

means <- colMeans(toxicity)
toxicity<- data.frame(scale(toxicity, center = T, scale = F))
x <- model.matrix(tox ~ C1 + C2 + C3 + C4 + C5 + C6)[,-1]
y <- tox

lambda.grid <- 10^seq(0,-2,length=100)
#lambda.grid <- seq(0.01,1,length.out = 100)

cv.lasso <- cv.glmnet(x,y,lambda=lambda.grid)
#cv.lasso <- cv.glmnet(x,y,lambda=lambda.grid,nfolds = length(y)) L1O

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

plot(cv.lasso)
abline(v=log(bestlam.lasso), lty=1)

fit.lasso.opt <- glmnet(x,y, lambda = bestlam.lasso) # default: alpha=1 -> lasso 
#we get the coefficients:
coef.lasso.opt <- predict(fit.lasso.opt, type = 'coefficients')
coef.lasso.opt
#So we select variables C1, C3, C4, C6

g1 <- lm(tox ~ C1 + C3 + C4 + C6)
summary(g1)

#verification of the assumptions:
shapiro.test(g$residuals)
x11()
plot(g$fitted.values, g$residuals/summary(g)$sigma) #few bad residuals, but still acceptable
anova(g,g1)


betas1 <- g1$coefficients
s <- summary(g1)$sigma

predict(g1, z_new, interval = "prediction", level = 1 - alpha)
