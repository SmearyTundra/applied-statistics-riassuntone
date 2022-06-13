setwd("Documents/Polimi/APPSTATS")

kimono <- read.table("kimono.txt")
# cast as factors
kimono$city <- as.factor(kimono$city)
kimono$type <- as.factor(kimono$type)
kimono$new = as.factor(paste(kimono$city, kimono$type))
L <- levels(kimono$new)

attach(kimono)

# verify normality of value for each subgroup
Pvals <- c(shapiro.test(value[ new==L[1] ])$p,
        shapiro.test(value[ new==L[2] ])$p,
        shapiro.test(value[ new==L[3] ])$p,
        shapiro.test(value[ new==L[4] ])$p)
Pvals # taking confidence level .95 it is passed

# homoschedasticity: visual inspection
x11()
boxplot(cbind(value[ new==L[1] ],value[ new==L[2] ],value[ new==L[3] ],
              value[ new==L[4] ]))
# we see the length of the whiskers of the boxes are the same practically.

# proceed with anova
fit <- aov(value  ~ city + type + city:type)

summary(fit)

# take off interaction as it is not significant
fit <- aov(value  ~ city + type)
summary(fit)

# now take also city as it is not significant either
fit <- aov(value  ~type)  # fit anova

L <- levels(as.factor(type))  # define levels
n1 <- length(value[type==L[1]])
n2 <-  length(value[type==L[2]])
mu.hat1 <- mean(value[type==L[1]])
mu.hat2 <- mean(value[type==L[2]])
var.1 <- var(value[type==L[1]])
var.2 <- var(value[type==L[2]])

# Perform confidence intervals
alpha = 0.05
# number of simultaneous tests: equal to the number of groups
# a1 will be [1, 0], a2 will be [0, 1] 
k = 2  
cfr.t <- qt(1-alpha/(2*k),n1+n2-2)
IC.mu.1 <- cbind(inf = mu.hat1 - cfr.t*sqrt(var.1/n1),
                 center =  mu.hat1,
                  sup = mu.hat1 + cfr.t*sqrt(var.1/n1))
IC.mu.2 <- cbind(inf = mu.hat2 - cfr.t*sqrt(var.2/n2),
                 center =  mu.hat2
                 , sup = mu.hat2 + cfr.t*sqrt(var.2/n2))
# TODO curiosità: taking linear combination of mu1 and mu2, by summing variances
# do we get the same as the pooled variance?

# varianza raggrupata: we know there is homoschedasticity so we take best estimate
var.raggr <- (var.1 * (n1-1) + var.2 * (n2-1)) / (n1+n2-1)
IC.diff <- cbind(inf = mu.hat2 - mu.hat1 - cfr.t*sqrt(var.raggr/(n1+n2-2)),
                 center =  mu.hat2 - mu.hat2,
                 sup =  mu.hat2 - mu.hat1 +  cfr.t*sqrt(var.raggr/(n1+n2-2)))



#############################33
### PROBLEM 2
###########################
bento <- read.table("bento.txt")
dim(bento)
n <- dim(bento)[1]
names(bento)
# get differences
D <- bento[,1:4] - bento[5:8]
# verify assumption of D


p <- dim(D)[2]

X.bar <- sapply(D, mean)
Sdiff <- cov(D)

delta.0 <- rep(0, p)
Tsq <- n * t((X.bar - delta.0)) %*% solve(Sdiff) %*%(X.bar - delta.0)
alpha = 0.05
cfr.F <-  (n-1)*p/(n-p) *qf(1-alpha,p,n-p)
Tsq < cfr.F

# p value: 
P <- 1-pf(Tsq * (n-p) / (n-1) / p,  p, n - p)
P

## 
# simultaneous confidence intervals
###

# with the F:                                                          
for(i in seq(1, p)){
  CI <- cbind( inf = X.bar[i] - sqrt( (n-1)*p / (n-p) * cfr.F )* sqrt(Sdiff[i,i]/n),
               med = X.bar[i],
               sup = X.bar[i] +  sqrt( (n-1)*p / (n-p) * cfr.F )* sqrt(Sdiff[i,i]/n))
  print(CI)}
# Bonferroni
k <- p
cfr.t <- qt(1-alpha/2/k, n-1)
for(i in seq(1, p)){
  CI <- cbind( inf = X.bar[i] - cfr.t * sqrt(Sdiff[i,i]/n),
               med = X.bar[i],
               sup = X.bar[i] +   cfr.t * sqrt(Sdiff[i,i]/n))
  print(CI)}
#####

#### 
# c: Clustering
####


##
# D: Linear model
####
garden <- read.table("garden.txt")
garden
attach(garden)

# A. estimate parameters
fm <- lm(extension ~ 1 + carps + maple + cherry + stones)
summary(fm)
fm$coefficients

# A. verify assumptions: 
shapiro.test(fm$residuals)
"Pass normality test"
x11()
par(mfrow=c(2,2))
plot(fm)
"There seem to be no patterns."

# B.i. Statistical test to check if dependence on maple or cherry:
library(car)
names(garden)
C <- rbind(c(0,0,1,0,0), c(0,0,0,1,0))
linearHypothesis(fm, C)
"We see it is significant"

# B.ii Statistical test to check if dependence of mean garden extension on lakes
C <- rbind(c(0,1,0,0,0), c(0,0,0,0,1))
# TODO: linearHypothesis(fm, C, c(0,0,0,0,0)) 
linearHypothesis(fm, C) 
"Significant"

# Comment on possible model weaknesses based on B.
"We see that it is very unlikely that either the lakes or the trees have no
effect, so we rejected both null hypotheses.
Given that when performing the overall regression, no coefficient is significant
taking into account we are making more than one test at a time, we try to 
shrink the model"
"The regressor with highest p value is discarded"
fm2 <- lm(extension ~ 1 + carps + maple + stones)
summary(fm2)
"We see now what the hypotheses on B assured us: it is very unlikely that neither
from the lake or tree variables are significant. By taking out the regressor
of cherry, we have now significant maple and stones as regressors "
"We proceed to shrink the model"
fm3 <- lm(extension ~ 1  + maple + stones)
summary(fm3)

fm3 <- lm(extension ~  maple + stones -1)
summary(fm3)
"We leave the model without intercept, although vector 1 does not belong 
any more to span of Z, and R squared loses interpretability."
"The updated estimated coefficients is 0 for all regressors except for
maple and stones, as seen in the last summary."

