leaven <- read.table("leaven.txt", header = T)
attach(leaven)
D <- ifelse(yeast=="sd",0,1)#D = 1 => by
B <- ifelse(yeast=="sd",1,0)#B = 1 => sd
#we write a model:
# volume = beta0 + {beta1_by*time or beta1_sd*time} + {beta2_by*time^2 or beta2_sd*time^2}

#volume = beta0 + beta1_by*time + beta2_by*time^2
#volume = beta0 + beta1_sd*time + beta2_sd*time^2

#volume = beta0 + beta1_g*time + beta2_g*time^2

g <- lm(volume ~ time + D:time + I(time^2) + D:I(time^2))
summary(g)
g1 <- lm(volume ~ time + B:time + I(time^2) + B:I(time^2))
summary(g1)

#verification of assumptions
shapiro.test(g$residuals)
x11()
qqnorm(g$residuals)
qqline(g$residuals)
x11()
plot(g$fitted.values, g$residuals/summary(g)$sigma)
abline(h = 2, col = 'red')
abline(h = -2, col = 'red')
#not ideal but the apparent heteroschedasticity might just caused by different number of observations across the fitted values
x11()
plot(time, g$residuals/summary(g)$sigma)
x11()
boxplot(g$residuals/summary(g)$sigma ~ yeast)
var.test(g$residuals ~ yeast)
#all assumptions are verified

betas <- c(beta0 = g$coefficients[1],
           beta1_by = g$coefficients[2] + g$coefficients[4],
           beta1_sd = g$coefficients[2],
           beta2_by = g$coefficients[3] + g$coefficients[5],
           beta2_sd = g$coefficients[3])
s <- summary(g)$sigma


betas1 <- c(beta0 = g1$coefficients[1],
           beta1_by = g1$coefficients[2],
           beta1_sd = g1$coefficients[2] + g1$coefficients[4],
           beta2_by = g1$coefficients[3],
           beta2_sd = g1$coefficients[3] + g1$coefficients[5])



an <- aov(volume ~ yeast)
summary(an)

#interpreting "rising of the dough" as the overall increase in volume, withouth accounting for the time spent
#Or we can interpret it also as:
#H0: beta1_sd - beta1_by = 0 & beta2_sd - beta2_by = 0 vs H1
#H0: beta_time:c = 0 & beta_time^2:c = 0
linearHypothesis(g, rbind(c(0,0,0,1,0), c(0,0,0,0,1)), c(0,0))
#indeed we have evidence to say that at least in one amongst the linear or quadratic behaviour, there is a significant difference between the two yeasts



#we want to simultaneously check:
#H0: beta2_sd = 0 & beta2_by =/=0 vs H1
#se we evaluate simultaneous confidence intervals, with bonferroni correction
library(car)
linearHypothesis(g, c(0,0,1,0,0),0)
linearHypothesis(g, c(0,0,1,0,1),0)
#assumptions for the test are already verified since we've proven that the model is gaussian

#we can reduce the model ignoring the quadratic term for sordough
#volume = beta0 + beta1_g*time + beta2_by*time^2
g_new <- lm(volume ~ time + D:time + D:I(time^2))
summary(g_new)
anova(g_new, g) #we are not losing information
shapiro.test(g_new$residuals)
x11()
plot(g_new$fitted.values, g_new$residuals/summary(g_new)$sigma)

betas_new <- c(beta0 = g_new$coefficients[1],
               beta1_by = g_new$coefficients[2] + g_new$coefficients[3],
               beta1_sd = g_new$coefficients[2],
               beta2_by = g_new$coefficients[4],
               beta2_sd = 0)
s_new <- summary(g_new)$sigma

#we generate confidence intervals (one at the time!) for both yeasts using this last model
alpha <- 0.01
z_sd <- data.frame(time = 2, D = 0)
z_by <- data.frame(time = 2, D = 1)
#we want to test:
#H0: beta1_by*2 + beta2_by*4 > beta1_sd*2 -> D:time + 2*D:time^2 > 0
linearHypothesis(g_new,c(0,0,1,2), mu = 0 )
#CI for D:time + 2*D:time^2 at 95%
a <- c(0,0,1,2)
cfr.stud <- qt(1-5*alpha, g_new$df)
Z <- cbind(rep(1, length(time)), time, D:time, D:I(time^2))
CI <- c(a%*%g_new$coefficients - cfr.stud*sqrt(s_new*a%*%solve(t(Z)%*%Z)%*%a),
        a%*%g_new$coefficients + cfr.stud*sqrt(s_new*a%*%solve(t(Z)%*%Z)%*%a))


predict(g_new, z_sd, interval = "confidence", level = 1-alpha)

#she will be better off using sourdough yeasts, since it yields better results with a confidence of 95%