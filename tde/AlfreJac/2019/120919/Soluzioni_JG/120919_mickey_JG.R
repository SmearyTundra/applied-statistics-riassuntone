mickey <- read.table("mickey.txt", header = T)
attach(mickey)
g <- lm(waiting.time ~ day.of.the.week + I(1 + cos(4*pi*day/365)) + day.of.the.week:I(1 + cos(4*pi*day/365)))
summary(g)
#we need to check gaussianity of the residuals
shapiro.test(g$residuals)
x11()
qqnorm(g$residuals)
qqline(g$residuals)
#we have that
x11()
plot(g$fitted.values,g$residuals/summary(g)$sigma)
x11()
boxplot(g$residuals ~ day.of.the.week)
x11()
plot(I(1 + cos(4*pi*day/365)), g$residuals/summary(g)$sigma)
#we do not see patterns
params <- c(alpha_0 = g$coefficients[1],
            alpha_1 = g$coefficients[1] + g$coefficients[2],
            beta_0 = g$coefficients[3],
            beta_1 = g$coefficients[3] + g$coefficients[4])
s <- summary(g)$sigma
params
#so on the weekend there is, on average, more waiting time, more or less with the same periodic behaviour of the (systemic) variability (1 + cos(4 * pi * day/365))
#Test:
#H0: alpha_1 - alpha_0 = 0 vs H1
library(car)
linearHypothesis(g, c(0,1,0,0), 0)
#we have strong evidence to reject the null hypothesis, p_val <<< 0.05
summary(g)
#by applying bonferroni correction  on the one-at the time test provided by summary(g), we can neglect at overall level 0.95 the impact of
#the interaction between categorycal variable and periodic one
#i.e. we have no evidence to say that beta_1 differs from beta_0
#hence we propose the reduced model:
# y = alpha_g + beta*(1 + cos(4 * pi * day/365)) + eps
g1 <- lm(waiting.time ~ day.of.the.week + I(1 + cos(4*pi*day/365)))

summary(g1)
anova(g,g1)
shapiro.test(g1$residuals)
x11()
plot(g1$fitted.values, g1$residuals)
params_1 <- c(alpha_0 = g1$coefficients[1],
                alpha_1 = g1$coefficients[1] + g1$coefficients[2],
                beta = g1$coefficients[3])
s_1 <- summary(g1)$sigma

#we have already strong evidence to support the fact that the maximum is achieved during weekends (one at the time confidence interval from summary)
#so we simultaneously (Bonferroni) test
linearHypothesis(g1, c(1,1,2), 60)
#we cannot reject the hypothesis:
#alpha_1 = 60 - 2*beta

#So we propose the reduced model
# y = alpha_0 + beta(1 + cos(4 * pi * day/365)) + I(g = 1)*(60 - 2*beta) + eps
W <- ifelse(day.of.the.week=="weekend",1,0)
wait <- waiting.time
wait[which(W==1)] <- wait[which(W==1)] - 60
g_red <- lm(wait ~ I(1-W) + I((1 + cos(4 * pi * day/365)) - 2*W) - 1)
summary(g_red)
shapiro.test(g_red$residuals)
params_red <- c(alpha_0 = g_red$coefficients[1],
                alpha_1 = 60 - 2*g_red$coefficients[2],
                beta = g_red$coefficients[2])
s_red <- summary(g_red)$sigma
z_new <- data.frame(W = 0, day = 238)
predict(g_red, z_new, interval = "prediction", level = 0.95)+60*(z_new$W==1)
#this result assumes gaussianity of the residuals (checked) and independence of the residuals from any regressor used
#i.e:
x11()
plot(day, g_red$residuals)
x11()
boxplot(g$residuals ~ day.of.the.week)
#homoschedasticity wrt to weekends vs weekdays is a bit week, but this could easily be imputed to the different sampling sizes
summary(as.factor(day.of.the.week))
