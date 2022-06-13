garden <- garden[,-1]#get rid of tag
attach(garden)
#create the model
n <- dim(garden)[1]
g <- lm(extension ~ carps + maple + cherry + stones)
summary(g)
betas <- g$coefficients
#verification of gaussianity of residuals
shapiro.test(g$residuals)
qqnorm(g$residuals)
qqline(g$residuals)
#gaussianity holds
#homoschedasticity of residuals and lack of patterns
x11()
plot(g$fitted.values,scale(g$residuals, scale = T, center = F))
points(g$fitted.values, rep(2,n), type = 'l', col = 'red')
points(g$fitted.values, rep(-2,n), type = 'l', col = 'red')
#all good, only a few important residuals, not a problem since we have a quite large n

x11()
plot(carps,scale(g$residuals, scale = T, center = F))
points(carps, rep(2,n), type = 'l', col = 'red')
points(carps, rep(-2,n), type = 'l', col = 'red')
#all good

x11()
plot(maple,scale(g$residuals, scale = T, center = F))
points(maple, rep(2,n), type = 'l', col = 'red')
points(maple, rep(-2,n), type = 'l', col = 'red')
#all good

x11()
plot(cherry,scale(g$residuals, scale = T, center = F))
points(cherry, rep(2,n), type = 'l', col = 'red')
points(cherry, rep(-2,n), type = 'l', col = 'red')
#all good

x11()
plot(stones,scale(g$residuals, scale = T, center = F))
points(stones, rep(2,n), type = 'l', col = 'red')
points(stones, rep(-2,n), type = 'l', col = 'red')
#all good

#residual assumptions hold


#we do NOT need to apply Bonferroni correction (or any SimCI correction), they would accounto for "they both have an impact, simultaneously", we just need that either one has an effect
#test_i : H0 : beta(i)= 0 vs H1
#T-statistic : abs(beta^(i)- 0)/(S*sqrt((Z'Z)^(-1)[i,i]) ~ t(n-(r+1))
#but here summary(g) already provides the results of this test (ONE AT THE TIME)
#
#The test required is og the kind:
#test_i : H0 : (beta(maple)= 0 & beta(cherry)= 0) vs H1
#we want to reject H0 if BOTH of the one at the time reject their individual H0

#maple -> 0.06 cherry -> 0.28   neither of them reject, but, keeping in mind that these p_values are quite conservative on the assumption,
#when we apply Bonferroni we divide the p_values for keeping into account that we want them to hold simultaneously
#the correction that would be needed here is not a simple multiplicative factor (P(union) = sumP(), while P(inters) = prodP(), but under independence, more complex computations would be needed)
#anyway the correction could only correct the p-values upward, so the fact that the p-value of cherry is border-line is not a problem

#carps->0.26 stones ->0.26, neither of them reject, trivially

#So we have a model which explains 65% of the variability, but which has no statistical evidence to say that any of the regressors used has an effect
#Not surprising, indeed computing the vif:
vif(g)
#we notice that we have a high collinearity of our regressors
#we can apply backward selection, excluding from the model the regressor with the highest p-value for the test of its coefficient being zero
#ONE AT THE TIME!
summary(g)
#we remove cherry
g1 <- lm(extension ~ carps + maple + stones)
summary(g1)
#we remove carps
g2 <- lm(extension ~ maple + stones)
summary(g2)
#now we have a model with a high significance of our regressors without losing in terms of p-value
#this model is valid, indeed
anova(g,g2)
#p_val = 0.5 -> the reduced model doesn't lose wrt the complete one
shapiro.test(g2$residuals)
x11()
plot(g2$fitted.values,scale(g2$residuals, scale = T, center = F))
points(g2$fitted.values, rep(2,n), type = 'l', col = 'red')
points(g2$fitted.values, rep(-2,n), type = 'l', col = 'red')
#all assumptions are still valid


#we notice that we still have in our model one variable about the lake elements (stones) and one about the trees (maple)
#there is a high pair-wise correlation of the regressors, which caused a high variability on our estimates and "masked" the actual impact of these regressors

betas_new <- c(g2$coefficients[1], 0, g2$coefficients[2], 0, g2$coefficients[3])
names(betas_new) <- c("b0", "b1", "b2", "b3", "b4")
betas_new