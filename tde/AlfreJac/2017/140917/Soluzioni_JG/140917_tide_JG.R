tide <- read.table("tide.txt", header = T)
t_0 <- 82
attach(tide)
n <- dim(tide)[1]
g <- lm(h ~ I(sin(2*pi*t/28)) + I(sin(pi*(t - t_0)/365)) + t)
summary(g)
coeff <- g$coefficients
sigma <- sqrt(sum((g$residuals)^2)/(g$df.residual))
summary(g)$sigma
#we check gaussianity of the residuals and homoschedasticity
shapiro.test(g$residuals)
plot(g$fitted.values,g$residuals/summary(g)$sigma)
#we test H0: {beta1 = 0 && beta2 =0} vs H1
#we apply a Bonferroni Correction, hence looking at the (one at the time) p-values of summary(g):
#p_val1 <- <2e-16, p_val2 <- 0.5130/2 = 0.2565
#we can say that there is dependence on periodic components, although to only one of them

#we test H0: beta0 = 0 vs H1
#we do not need to apply a correction <- p_value < 2e-16
#there is dependence

g1 <- lm(h ~ I(sin(2*pi*t/28)) + t)
summary(g1)

shapiro.test(g1$residuals)
sigma1 <- sqrt(sum((g1$residuals)^2)/(g1$df))
coeff1 <- g1$coefficients

z_new <- data.frame(t = c(263, 335))
#we probably want both prediction to hold simultaneously, we need to adjust alpha = 0.1/2 = 0.05

h_pred <- predict(g1, z_new, interval='prediction', level = 0.95)
rownames(h_pred)<- z_new$t
h_pred
#we can say with overall confidence 90% that on 1st Dec we won't have high water in Venice, while we can expect it on Sep the 20th