airport <- read.table("airport.txt", header = T)
attach(airport)
time <- as.factor(time.of.the.day)
g <- lm(duration ~ time + distance + time:distance)

summary(g)
#verification of the assumptions
shapiro.test(g$residuals)
x11()
qqnorm(g$residuals)
qqline(g$residuals)
#we can assume gaussianity
x11()
plot(g$fitted.values, g$residuals/summary(g)$sigma)
x11()
boxplot(g$residuals ~ time)
bartlett.test(g$residuals ~ time)
#so we have homoschedasticity amongst groups
x11()
plot(distance, g$residuals/summary(g)$sigma)
#the assumptions are verified

params <- c(beta0_6 = g$coefficients[1] + g$coefficients[3],
            beta0_11 = g$coefficients[1],
            beta0_16 = g$coefficients[1] + g$coefficients[2],
            beta1_6 = g$coefficients[4] + g$coefficients[6],
            beta1_11 = g$coefficients[4],
            beta1_16 = g$coefficients[4] + g$coefficients[5])
params
s <- summary(g)$sigma
s
#so we want to test:
#H0: beta0_6 = beta0_11 = beta0_16 & beta1_6 = beta1_11 = beta1_16 vs H1
library(car)
linearHypothesis(g, rbind(c(0,1,0,0,0,0),c(0,0,1,0,0,0),c(0,0,0,0,1,0),c(0,0,0,0,0,1)), c(0,0,0,0))
#There is strong evidence to reject the null hypothesis, there is a significant difference between at least two groups in at least one coefficient

#H0: beta1_6 = beta1_11 = beta1_16 = 0 vs H1
linearHypothesis(g, rbind(c(0,0,0,1,0,0),c(0,0,0,1,1,0),c(0,0,0,1,0,1)), c(0,0,0))
#there is strong evidence to reject the null hypothesis, in at least one group there is a significant dependence on the distance travelled

summary(g)

#we can start by tring to neglect the dependence on time only in the intecept
#H0: beta0_6 = beta0_11 = beta0_16 vs H1
linearHypothesis(g, rbind(c(0,1,0,0,0,0),c(0,0,1,0,0,0)), c(0,0))
#there is no evidence to rejct the null hypothesis

#we start proposing a reduced model with dependence on time of the day only on the distance travelled
g1 <- lm(duration ~ distance + time:distance)
summary(g1)

anova(g,g1)
shapiro.test(g1$residuals)
x11()
plot(g1$fitted.values, g1$residuals)
#still valid and still significant

g$coefficients
#we can try and test if there is no difference on the coefficients for distance between 16-20 and 11-15
#H0: beta1_16 = beta1_11 vs H1
#and by looking at the (one at the time) p_value for this test on the summary we can say that there is not

#seeing if we can neglect all 4 of them from the general model:
linearHypothesis(g, rbind(c(0,1,0,0,0,0),c(0,0,1,0,0,0),c(0,0,0,0,1,0)), c(0,0,0))
#still no evidence to reject the null hypothesis

#Creating the new dummy variable
TM <- ifelse(time=="6-10", 1, 0) # 1 if in 6-10
g_red <- lm(duration ~ distance + TM:distance)

shapiro.test(g_red$residuals)
x11()
plot(g_red$fitted.values, g_red$residuals/summary(g)$sigma)
#still valid

params_red <- c(beta0 = g_red$coefficients[1],
                beta1_6 = g_red$coefficients[2] + g_red$coefficients[3],
                beta1_11 = g_red$coefficients[2],
                beta1_16 = g_red$coefficients[2])
s_red <- summary(g_red)$sigma

#using this last reduced model
z_new <- data.frame(TM = 1, distance = 57)
#we want a prediction interval for the duration of the trip, then using the upper estimate to take the best bus
predict(g_red, z_new, interval = "prediction", level = 0.99)
#so it takes at most 1 hour and a half, hence in order to be there at 9:30 he should take the bus leaving at 7 am,
#having still a 3 minute margin on the worst case scenario