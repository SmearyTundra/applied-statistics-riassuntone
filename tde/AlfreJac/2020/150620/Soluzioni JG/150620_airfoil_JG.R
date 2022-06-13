airfoil <- read.table("airfoil.txt", header = T)
attach(airfoil)
#create dummy
H <- ifelse(velocity == "H",1,0)#H = 1 => high velocity
g <- lm(sound ~ H + frequency + H:frequency)
summary(g)
shapiro.test(g$residuals)
x11()
qqnorm(g$residuals)
qqline(g$residuals)
#we can work with gaussianity
x11()
plot(g$fitted.values, g$residuals/summary(g)$sigma)
x11()
boxplot(g$residuals ~ velocity)
x11()
plot(frequency, g$residuals/summary(g)$sigma)
#we can assume homoschedasticity

betas <- c(beta0_H = g$coefficients[1] + g$coefficients[2],
           beta0_L = g$coefficients[1],
           beta1_H = g$coefficients[3] + g$coefficients[4],
           beta1_L = g$coefficients[3])
s <- summary(g)$sigma


#we need to perform the TESTs
#H0: beta1_H = 0 & beta1_L = 0  vs H1    so either case the mean sound level is not influenced by frequency
library(car)
linearHypothesis(g, rbind(c(0,0,1,0), c(0,0,1,1)), c(0,0))
#There is a significan impact in at least one of the two scenarios

#H0: beta0_H = 0 & beta0_L = 0  vs H1    so either case the mean sound level is not influenced by air velocity
linearHypothesis(g, rbind(c(1,0,0,0), c(1,1,0,0)), c(0,0))
#There is a significan impact in at least one of the two scenarios


#H0: beta1_H = beta1_L = 0  vs H1
linearHypothesis(g, c(0,0,0,1),0)
#we do not have evidence to reject the null hypothesis

#so we propose the reduced model:
# Y = beta0_g + beta1*x + eps
g_new <- lm(sound ~ H + frequency)
summary(g_new)

anova(g,g_new)
shapiro.test(g_new$residuals)
x11()
plot(g_new$fitted.values, g_new$residuals/summary(g_new)$sigma)
x11()
boxplot(g_new$residuals~velocity)
x11()
plot(frequency, g_new$residuals/summary(g_new)$sigma)
#all assumptions are still verified and we do not lose much considering just the reduced model
betas_new <- c(beta0_H = g_new$coefficients[1] + g$coefficients[2],
               beta0_L = g$coefficients[1],
               beta1 = g$coefficients[3])
s_new <- summary(g_new)$sigma

#it is a confidence interval
z_new <- data.frame(frequency = 15000, H = 1)
alpha <- 0.05
predict(g_new, z_new, level = 1 - alpha, interval = "confidence")
