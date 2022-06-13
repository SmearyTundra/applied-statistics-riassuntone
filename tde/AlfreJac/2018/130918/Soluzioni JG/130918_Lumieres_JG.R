Lumieres <- read.table("Lumieres.txt", header = T)
attach(Lumieres)
R <- ifelse(rain == "yes",1, 0)
g <- lm(N.people ~ R + day + I(day^2) + temperature)

summary(g)
betas <- c(beta0_1 = g$coefficients[1],
           beta0_1 = g$coefficients[1] + g$coefficients[2],
           beta1 = g$coefficients[3],
           beta2 = g$coefficients[4],
           beta3 = g$coefficients[5])
s <- summary(g)$sigma
#gaussianity:
shapiro.test(g$residuals)
x11()
qqnorm(g$residuals)
qqline(g$residuals)
#gaussianity ok
x11()
plot(g$fitted.values, g$residuals/s)
x11()
plot(day, g$residuals/s)
x11()
plot(I(day^2), g$residuals/s)
x11()
plot(temperature, g$residuals/s)
#everything is ok
#no suspects of leverages

summary(g)
#given the one at the time p_values for the tests:
#H0: beta1 = 0 vs H1
#H0: beta2 = 0 vs H1
#we can see that there is strong evidence supporting a dependence of the number of partecipant on the day, both linearly and quadratically


#we perform a test on the coefficients of beta, simultaneously checking:
#H0: beta0_1 - beta0_2 = 0 & beta3 = 0
library(car)
linearHypothesis(g, rbind(c(0,1,0,0,0), c(0,0,0,0,1)), c(0,0))
#the p_value is much higher than 0.5, so we cannot reject the null hypothesis and we propose the reduced model:
#n = beta0 + beta1*d + beta2*d^2 + eps
g_new <- lm(N.people ~ day + I(day^2))



anova(g, g_new)
shapiro.test(g_new$residuals)
x11()
plot(g_new$fitted.values, g_new$residuals/summary(g_new)$sigma)


betas_new <- c(beta0 = g_new$coefficients[1],
               beta1 = g_new$coefficients[2],
               beta2 = g_new$coefficients[3])
s_new <- summary(g_new)$sigma

#we are testing:
#H0: beta_1 + 2*beta_2*61 = 0 vs H1

linearHypothesis(g_new, c(0,1,2*61), 0)
#we do not have evidence to reject H0
#beta_1 = -122beta_2
#So we can further reduce the model:
#n = beta0 +beta(-122*d + d^2) + eps


#n = beta_0 + beta*(d^2 - 122*d) + eps
g_red <- lm(N.people ~ I(day^2 - 122*day))

anova(g_red,g)
shapiro.test(g_red$residuals)
x11()
plot(g_red$fitted.values, g_red$residuals/summary(g_red)$sigma)

betas_red <- c(beta0 = g_red$coefficients[1],
               beta = g_red$coefficients[2])
s_red <- summary(g_red)$sigma

summary(g_red)


z_new <- data.frame(day = 58, temperature = 29, R = 0)
alpha <- 0.05
predict(g_red, z_new, interval = "prediction", level = 1- alpha)
