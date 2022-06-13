albatross <- read.table("albatross.txt", header = T)
n <- dim(albatross)[1]
#dummy -> 1 if upwind
D <- rep(0, n)
D[which(albatross$wind == "upwind")] <- rep(1, length(which(albatross$wind == "upwind")))
attach(albatross)
g <- lm(distance ~ D + I(Va^2) + I(Vi^2) + D:I(Va^2) +  D:I(Vi^2))
summary(g)
#alpha1 = Int + D
al1 <- g$coefficients[1] + g$coefficients[2]
#alpha2 = Int
al2 <- g$coefficients[1]
#beta1 = D:Va^2 + Va^2
be1 <- g$coefficients[3] + g$coefficients[5]
#beta2 = Va^2
be2 <- g$coefficients[3]
#gamma1 = D:Vi^2 + Vi^2
ga1 <- g$coefficients[4] + g$coefficients[6]
#gamma2 = Vi^2
ga2 <- g$coefficients[4]
si <- summary(g)$sigma

shapiro.test(g$residuals)
plot(g$fitted.values,g$residuals/si)

#we reduce by eliminating the non significant regressors -> one at the time!
summary(g)
g1 <- lm(distance ~ I(Va^2) + I(Vi^2) + D:I(Va^2) +  D:I(Vi^2)) #same average
summary(g1)
anova(g1, g) #we do not lose with the reduced model
shapiro.test(g1$residuals)
plot(g1$fitted.values, g1$residuals/summary(g1)$sigma)

#simultaneous test:
#H0 : {gamma1 = - beta1 && gamma2 = - beta2} vs H1
library(car)
linearHypothesis(g1, rbind(c(0,1,1,0,0), c(0,1,1,1,1)), c(0,0))
#or with Bonferroni
linearHypothesis(g1, rbind(c(0,1,1,0,0)), 0)
linearHypothesis(g1, rbind(c(0,1,1,1,1)), 0)
#anyway we cannot reject the null hypothesis
#we can use the constrained model:
g2 <- lm(distance ~ I(Va^2 - Vi^2) + D:I(Va^2 - Vi^2)) #since the coefficients of Vi^2 would be - coeff of Va^2
summary(g2)
anova(g1, g2) #we do not lose with the constrained model
alpha2_new <- alpha1_new <- g2$coefficients[1]
beta1_new <- g2$coefficients[2] + g2$coefficients[3]
gamma1_new <- - beta1_new
beta2_new <- g2$coefficients[2]
gamma2_new <- - beta2_new

#we need Bonferroni correction, they are PREDICTION intervals:
alpha <- 0.01/2
z_upwind <- data.frame(Va = 35, Vi = 25, D = 1)
z_downwind <- data.frame(Va = 35, Vi = 25, D = 0)
predict(g2, z_upwind, interval = 'prediction',level = 1-alpha)
predict(g2, z_downwind, interval = 'prediction',level = 1-alpha)
#upwind is safe, downwind is not