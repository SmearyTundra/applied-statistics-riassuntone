kimonos <- read.table("kimonos.txt", header = T)
attach(kimonos)

g <- lm(price ~ q + l + q:l + ncolors + q:ncolors)

shapiro.test(g$residuals)
x11()
plot(g$fitted.values, g$residuals/summary(g)$sigma)
x11()
boxplot(g$residuals/summary(g)$sigma ~ H)
#all good
params <- c(alpha_2 = g$coefficients[1] + g$coefficients[2],
            alpha_1 = g$coefficients[1],
            beta_2 = g$coefficients[3] + g$coefficients[5],
            beta_1 = g$coefficients[3],
            gamma_2 = g$coefficients[4] + g$coefficients[6],
            gamma_1 = g$coefficients[4])
params
s <- summary(g)$sigma

summary(g)

library(car)
#test:
#H0: alpha_1 = alpha_2 vs H0
#as reported in the summary the p_val for the one at the time test for param_H = 0 is 0.62
#we cannot reject the null hypothesis

#test:
#H0: gamma_1 = 0 & gamma_2 = 0 vs H1
linearHypothesis(g, rbind(c(0,0,0,1,0,0),c(0,0,0,1,0,1)), c(0,0))
#we can reject at level 10%

g1 <- lm(price ~ l + q:l + ncolors + q:ncolors)

shapiro.test(g1$residuals)
x11()
plot(g1$fitted.values, g1$residuals/summary(g1)$sigma)

params1 <- c(alpha = g1$coefficients[1],
            beta_2 = g1$coefficients[2] + g1$coefficients[4],
            beta_1 = g1$coefficients[2],
            gamma_2 = g1$coefficients[3] + g1$coefficients[5],
            gamma_1 = g1$coefficients[3])
params1
s1 <- summary(g1)$sigma

summary(g1)

#one step more, in one of the two groups n.colors is descardable
#which one?

g2 <- lm(price ~ l + q:l + q:ncolors)
summary(g2)

shapiro.test(g2$residuals)
x11()
plot(g2$fitted.values, g2$residuals/summary(g2)$sigma)

summary(g2)

#dummy variable
M <- ifelse(q=="medium",1,0)
g3 <- lm(price ~ l + M:l + M:ncolors)
summary(g3)

shapiro.test(g3$residuals)
x11()
plot(g3$fitted.values, g3$residuals/summary(g3)$sigma)


params_red <- c(alpha = g3$coefficients[1],
                beta_1 = g3$coefficients[2],
                beta_2 = g3$coefficients[2] + g3$coefficients[3],
                gamma_2 = g3$coefficients[4],
                gamma_1 = 0)
params_red
s_red <- summary(g3)$sigma


#since we neglected the effect on the mean of the quality, we are activly looking for 3 confidence intervals, Bonferroni correction:
alphaB <- 0.1/3
alpha <- data.frame(l=0, M = 0, ncolors = 0)
z_new_m <- data.frame(l=6.5, M = 1, ncolors = 3)
z_new_h <- data.frame(l=6.5, M = 0, ncolors = 3)

predict(g3, alpha,level = 1 - alphaB, interval = "confidence")
predict(g3, z_new_h,level = 1 - alphaB, interval = "confidence")
predict(g3, z_new_m,level = 1 - alphaB, interval = "confidence")
