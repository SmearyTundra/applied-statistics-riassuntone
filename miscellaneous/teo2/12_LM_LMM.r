# Topics:
#   - Linear Models with heteroschedastic/dependent errors
#   - Linear Mixed Models with Random Intercept
#   - Linear Mixed Models with Random Intercept + Random Slope
#   - Inference
#   - Assessing Assumptions
#   - Prediction



rm(list=ls())
graphics.off()


library(ggplot2)
library(insight)
library(lattice)
library(lme4)
library(nlme)
library(corrplot)
library(plot.matrix)

############# BEFORE STARTING ##############
# Importing dataset
school = read.table('school.txt', header=T)
school$gender= as.factor(school$gender)
school$school_id= as.factor(school$school_id)
head(school)

# Visual exploration to see if it is worth using LMM
x11()
ggplot(data=school, aes(x=as.factor(school_id), y=achiev, fill=as.factor(school_id))) +
  geom_boxplot() +
  labs(x='Primary School', y='Achievement') +
  ggtitle('Boxplot of achievements among primary schools') +
  theme_minimal() +
  theme(axis.text=element_text(size=rel(1.15)),axis.title=element_text(size=rel(1.5)),
        plot.title = element_text(face="bold", size=rel(1.75)), legend.text = element_text(size=rel(1.15)),
        legend.position = 'none')


############# LINEAR MODELS ##############
# MODEL: achiev_i = beta_0 + beta_1*gender_i+ beta_2*escs_i + eps_i
# eps_i ~ N(0, sigma2_eps)

lm1 = lm(achiev ~ gender + escs + school_id, data = school)
summary(lm1)
plot(lm1$residuals)
boxplot(lm1$residuals ~ school$school_id, col='orange', xlab='School ID', ylab='Residuals')

# Boxplot residuals according to different group (useful for LMM later)
colori = rainbow(length(unique(armd$subject)))
boxplot(lm1.form$residuals ~ armd$subject, col=colori,
        xlab='Subjects', ylab='Residuals', main ='Distribution of residuals across patients')  ## --> informative!

# Color or boxplot residuals according to different quantity that may affect variance (useful for extended LM now)
set.seed(1)
colori = rainbow(4)
colori2 = colori[armd$tp] # associate to each one of the 4 time instants a color
plot(lm1.form$residuals, col=colori2, ylab='residuals')
abline(h=0)

boxplot(lm1.form$residuals ~ armd$time.f, col=colori,
        xlab='Time.f', ylab='Residuals')  ## -> the variance of th observations increases in time




##### with heteroscedastic and independent errors
## OPTION 1 : VarIdent()
fm9.1 <- gls(visual ~ -1 + visual0 + time.f + treat.f:time.f,  # the same as before
             weights = varIdent(form = ~1|time.f), # Var. function; <delta, stratum>-group
             data = armd)
summary(fm9.1)
plot(fm9.1$residuals) 
fm9.1$modelStruct$varStruct
intervals(fm9.1, which = "var-cov")  ## 95% CI
anova(fm9.1, lm1.form) 

## OPTION 2 : VarPower()
fm9.2 <- update(fm9.1, weights = varPower(form = ~time)) # Var. function; <delta, v_it>-group
summary(fm9.2)
fm9.2$modelStruct$varStruct
intervals(fm9.2, which = "var-cov")
anova(fm9.2, fm9.1)


##### with heteroscedastic and dependent errors 
## (heteroschedasticity built upon previous section)

## choose dependence structure using variogram

## Variogram per group (se ho un factor essenziale as.numeric)
Vg2 <- Variogram(fm9.2, form = ~ tp | subject)
Vg2
plot(Vg2, smooth = FALSE, xlab = "Time Lag",ylim=c(0,0.7))
# Il variogram aumenta -> la correlazione diminuisce -> structure AR(1)
# Il variogram non ha pattern -> compound symmetry -> CorCompSym()

## OPTION 1: CorCompSym()
lm1.form <- formula(visual ~ -1 + visual0 + time.f + treat.f:time.f )
fm12.1 <- gls(lm1.form, weights = varPower(form = ~time),
        correlation = corCompSymm(form = ~1|subject),
        data = armd)
summary(fm12.1)
intervals(fm12.1, which = "var-cov")

# Variance-Covariance structure
fm12.1vcov <- getVarCov(fm12.1, individual = "2")  #estimate of R_i, e.g. i=2
nms <- c("4wks", "12wks", "24wks", "52wks")
dnms <- list(nms, nms) # Dimnames created
dimnames(fm12.1vcov) <- dnms # Dimnames assigned
print(fm12.1vcov)
print(cov2cor(fm12.1vcov), corr = TRUE, stdevs = FALSE)  # Estimate of C_i (correlation matrix)

## OPTION 2: AR(1)
fm12.2 <- update(fm9.2, 
                 correlation = corAR1(form = ~tp|subject),
                 data = armd)

## OPTION 3: general correlation structure
fm12.3 <- update(fm12.2, correlation = corSymm(form = ~tp|subject),  ## the variance function is still VarPower()
                 data = armd)

## Final model fit:
panel.bwxplot0 <- function(x,y, subscripts, ...){
                        panel.grid(h = -1)
                        panel.stripplot(x, y, col = "grey", ...)
                        panel.bwplot(x, y, pch = "|", ...)
                        }
bwplot(resid(fm12.3) ~ time.f | treat.f, 
         panel = panel.bwxplot0,
         ylab = "Residuals", data = armd)
# Is the modeled heteroschedasticity present?

plot(fm12.3) 
# Is the modeled dependence present? Ex. highest residuals in some groups
# try plotting residuals for different groups and it shouldn't be present anymore

plot(fm12.3, 
       resid(., type = "p") ~ fitted(.) | time.f)










############# LINEAR MIXED EFFECTS MODELS (with intercept) ##############
# MODEL: achiev_ij = beta_0 + beta_1*gender_ij + beta_2*escs_ij + b_0i + eps_ij
# eps_ij ~ N(0, sigma2_eps)
# b_i ~ N(0, sigma2_b)
# random intercept -> school effect

lmm1 = lmer(achiev ~ gender + escs + (1|school_id), 
                      data = school)
summary(lmm1)

# Fixed Effects and 95% CIs
#-------------------------------
confint(lmm1, oldNames=TRUE)
fixef(lmm1)

# Variance components
#--------------------
print(vc <- VarCorr(lmm1), comp = c("Variance", "Std.Dev."))
sigma2_eps <- as.numeric(get_variance_residual(lmm1))
sigma2_eps
sigma2_b <- as.numeric(get_variance_random(lmm1))
sigma2_b

PVRE <- sigma2_b/(sigma2_b+sigma2_eps)
PVRE # intraclass correlation >20% high

# Random effects: b_0i
#----------------------------
ranef(lmm1)

# The dotplot shows the point and interval estimates for the random effects, 
# ordering them and highlighting which are significantly different from the mean (0)
# More informative groups -> more observation in this group -> smaller interval
x11()
dotplot(ranef(lmm1))

# Random intercepts and fixed slopes: (beta_0+b_0i, beta_1, beta_2)
coef(lmm1)
head(coef(lmm1)$school_id)

# Diagnostic plots 
#------------------
# 1) Assessing Assumption on the within-group errors
x11()
plot(lmm1)

x11()
qqnorm(resid(lmm1))
qqline(resid(lmm1), col='red', lwd=2)

# 2) Assessing Assumption on the Random Effects
x11()
qqnorm(unlist(ranef(lmm1)$school_id), main='Normal Q-Q Plot - Random Effects for Primary School')
qqline(unlist(ranef(lmm1)$school_id), col='red', lwd=2)

# Prediction
#-------------
# Prediction from regression model
predict_lm <- predict(lm1)
head(predict_lm)

# Prediction from mixed model:
# 1) Without random effects ->  re.form=NA
predict_no_re <- predict(lmm1, re.form=NA)
head(predict_no_re) # same predictions
# 2) With random effects
predict_re <- predict(lmm1)     ## --> remember to allow new levels in the RE if any
head(predict_re)

## Scenario Analysis
new_student1 = data.frame(gender=as.factor(1), escs=0.7, school_id=as.factor(32)) # observed school
new_student2 = data.frame(gender=as.factor(1), escs=0.7, school_id=as.factor(11)) # observed school
new_student3 = data.frame(gender=as.factor(1), escs=0.7, school_id=53) # new school

predict(lmm1, new_student1, re.form=NA)
predict(lmm1, new_student1)

predict(lmm1, new_student2, re.form=NA)
predict(lmm1, new_student2)

predict(lmm1, new_student3, re.form=NA)
predict(lmm1, new_student3, allow.new.levels = T) # student going in new school (level not observed)
# The prediction for new schools has to be the same with or without random effects

















############ LINEAR MIXED EFFECTS MODELS (with random intercept and slope) ##############
# Association between y and another regressor varies across groups

# MODEL:  achiev_ij = beta_0 + b_0i + (beta_1 + b_1i)*escs_i + eps_i --> homoscedastic residuals 
# random slope: b_1i*escs_i

# eps_i ~ N(0, sigma2_eps)
# Random effects: b_i ~ N(0, Sigma)

# To allow both the intercept, represented by 1, and the slope, represented by escs,
# to vary by student we can add the term:
#   - (1+escs|school_id)
# or, in alternative, without 1
#   - (escs|school_id)

lmm2 = lmer(achiev ~ gender + escs + (1 + escs|school_id), 
                data = school)
summary(lmm2)

# Fixed Effects and 95% CIs
#-------------------------------
confint(lmm2, oldNames=TRUE)
# Fixed effects: (beta_0, beta_1, beta_2)
fixef(lmm2)


# Variance components
#--------------------
print(vc <- VarCorr(lmm2), comp = c("Variance", "Std.Dev."))
sigma2_eps <- as.numeric(get_variance_residual(lmm2))
sigma2_eps
sigma2_b <- as.numeric(get_variance_random(lmm2))  ## it automatically computes Var(b0,b1)
sigma2_b


PVRE <- sigma2_b/(sigma2_b+sigma2_eps)
PVRE # intraclass correlation >20% high

# Random effects: b_0i
#----------------------------
ranef(lmm2) #(b_0i, b_1i) for i=1,...,200
head(ranef(lmm2)$school_id)

x11()
dotplot(ranef(lmm2))

# Random intercepts and fixed slopes: (beta_0+b_0i, beta_1, beta_2+b_2i)
coef(lmm2)
head(coef(lmm2)$school_id)


# Diagnostic plots 
#--------------------
# 1) Assessing Assumption on the within-group errors
x11()
plot(lmm2)

x11()
qqnorm(resid(lmm2))
qqline(resid(lmm2), col='red', lwd=2)


# 2) Assessing Assumption on the Random Effects
x11()
par(mfrow=c(1,2))
qqnorm(unlist(ranef(lmm2)$school_id[1]), main='Normal Q-Q Plot - Random Effects on Intercept')
qqline(unlist(ranef(lmm2)$school_id[1]), col='red', lwd=2)
qqnorm(unlist(ranef(lmm2)$school_id[2]), main='Normal Q-Q Plot - Random Effects on escs')
qqline(unlist(ranef(lmm2)$school_id[2]), col='red', lwd=2)

x11()
plot(ranef(lmm2))

anova(lmm1,lmm2)





########### VISUALIZATION OF RANDOM EFFECTS
## Visualization of random effects 
x11()
par(mfrow=c(1,3))
plot(c(1:50), unlist(coef(lmm2)$school_id[1]),
     xlab='School i', ylab=expression(beta[0]+b['0i']),
     pch=19, lwd=2, col='darkblue',
     main='Estimated random intercepts')
abline(h=fixef(lmm2)[1], lty=2, col='red', lwd=2)
legend(30, 13.5, legend=expression(paste('Fixed intercept ',beta[0])), lwd=2, lty=2, col='red', x.intersp=0.5)

plot(c(1:50), unlist(coef(lmm2)$school_id[2]),
     xlab='School i', ylab=expression(beta[1]),
     pch=19, lwd=2, col='darkblue',
     main='Estimated fixed slope for gender')
abline(h=fixef(lmm2)[2], lty=2, col='red', lwd=2)
legend(30,-0.6, legend=expression(paste('Fixed slope ',beta[1])), lwd=2, lty=2, col='red', x.intersp=0.5)

plot(c(1:50), unlist(coef(lmm2)$school_id[3]),
     xlab='Student i', ylab=expression(beta[2]+b['1i']),
     pch=19, lwd=2, col='darkblue',
     main='Estimated random slopes for escs')
abline(h=fixef(lmm2)[3], lty=2, col='red', lwd=2)
legend(30, 5, legend=expression(paste('Fixed slope ',beta[2])), lwd=2, lty=2, col='red', x.intersp=0.5)


# Lines Visualization
#---------------------
# Let's plot all the regression lines in 2D: if I have another categorical variable I do separate plots
## FEMALES
x11()
par(mfrow=c(1,2))
plot(school$escs[school$gender==0], school$achiev[school$gender==0],col='blue',
     xlab='escs', ylab='achievement',ylim=c(-5,30),main='Data and regression lines for females')
abline(10.0546535,1.6790886, col='red', lw=6)          

for(i in 1:50){
  abline(coef(lmm2)$school_id[i,1], coef(lmm2)$school_id[i,3])
}

## MALES
plot(school$escs[school$gender==1], school$achiev[school$gender==1],col='blue',
     xlab='escs', ylab='achievement',ylim=c(-5,30),main='Data and regression lines for males')
abline(10.02507-0.91180,1.96618, col='red', lw=6)  

for(i in 1:50){
  abline(coef(lmm2)$school_id[i,1] + coef(lmm2)$school_id[i,2], coef(lmm2)$school_id[i,3])
}

