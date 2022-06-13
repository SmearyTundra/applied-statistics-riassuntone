kimono <- read.table("kimono.txt",header = T)
attach(kimono)
#generate sub-datasets
kim_K_h <- kimono[which(city == "Kyoto" & type == "hand-made"),]
kim_T_h <- kimono[which(city == "Tokyo" & type == "hand-made"),]
kim_K_r <- kimono[which(city == "Kyoto" & type == "ready-to-use"),]
kim_T_r <- kimono[which(city == "Tokyo" & type == "ready-to-use"),]

av <- aov(value ~ city + type)

#gaussianity of residuals:
shapiro.test(av$residuals)
shapiro.test(av$residuals[which(city == "Kyoto" & type == "hand-made")])
shapiro.test(av$residuals[which(city == "Tokyo" & type == "hand-made")])
shapiro.test(av$residuals[which(city == "Kyoto" & type == "ready-to-use")])
shapiro.test(av$residuals[which(city == "Tokyo" & type == "ready-to-use")])
#homoschedasticity:
var.test(kim_K_h$value,kim_K_r$value)
var.test(kim_K_h$value,kim_T_h$value)
var.test(kim_K_h$value,kim_T_r$value)
var.test(kim_T_h$value,kim_T_r$value)
var.test(kim_T_h$value,kim_K_r$value)
var.test(kim_T_r$value,kim_K_r$value)


boxplot(value ~ city + type, data = kimono)

summary.aov(av)

#it suffices to see the P-value of the F test already shown in the summary to notice that the treatment of the city doesn't have an effect

av_new <- aov(value ~ type)
summary(av_new)

#verify that we can just use the reduced model
anova(av, av_new)

#we've reduced the model, there are now only two groups -> there is no need for a Bonferroni correction
#anyway a Bonferroni correction would have meant to divide the value of alpha by the number of simultaneous confidence intervals that we want to get (now just one!)
n <- length(kimono$value)/2
g <- 2
SSres <- sum(residuals(av_new)^2)
Means <- tapply(value, type, mean)
#Pivotal quantity -> delta_sample_mean/sqrt(2*S2/n) ~ t(n-1)
IC <- c(diff(Means) - qt(0.95, n*g-g) * sqrt(2*SSres/(n*(g*n-g))), 
        diff(Means) + qt(0.95, n*g-g) * sqrt(2*SSres/(n*(g*n-g))))
names(IC) <- c('Inf', 'Sup')
IC
