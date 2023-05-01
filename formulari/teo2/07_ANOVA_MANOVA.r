### TOPICS:
### One-way ANOVA and MANOVA
### Two-way ANOVA and MANOVA

## bartlett.test(d.numeric, as.factor(d$type1):as.factor(d$type2))

################# One-way ANOVA ##############

### Model: weigth.ij = mu + tau.i + eps.ij; eps.ij~N(0,sigma^2)
### Test:
### H0: tau.1 = tau.2 = tau.3 = tau.4 = tau.5 = tau.6 = 0
### H1: (H0)^c
### i.e.,
### H0: The feed supplements don't have effect
###     (= "chickens belong to a single population")
### H1: At least one feed supplement has effect
###     (= "chickens belong to 2, 3, 4, 5 or 6 populations")

# feed <- categorical, weight <- y

# Boxplot
plot(d$feed, d$weight, xlab='treat', ylab='weight', col='grey85', main='Dataset Chicken Weights')

# PLOT OF THE TWO MODELS
par(mfrow=c(1,2))
barplot(rep(mean(weight),6), names.arg=levels(feed), ylim=c(0,max(weight)),
        las=2, col='grey85', main='Model under H0')
barplot(tapply(weight, feed, mean), names.arg=levels(feed), ylim=c(0,max(weight)),
        las=2, col=rainbow(6),main='Model under H1')


n       <- length(d)      # total number of obs.
ng      <- table(d$feed)       # number of obs. in each group
treat   <- levels(d$feed)      # levels of the treatment
g       <- length(d$treat)     # number of levels (i.e., of groups)


#### VERIFYING ASSUMPTIONS
###-----------------
# 1) normality (univariate) in each group (6 tests)
Ps <- c(shapiro.test(weight[ feed==treat[1] ])$p,
        shapiro.test(weight[ feed==treat[2] ])$p,
        shapiro.test(weight[ feed==treat[3] ])$p,
        shapiro.test(weight[ feed==treat[4] ])$p,
        shapiro.test(weight[ feed==treat[5] ])$p,
        shapiro.test(weight[ feed==treat[6] ])$p) 
Ps

# 2) same covariance structure (= same sigma^2)
Var <- c(var(weight[ feed==treat[1] ]),
         var(weight[ feed==treat[2] ]),
         var(weight[ feed==treat[3] ]),
         var(weight[ feed==treat[4] ]),
         var(weight[ feed==treat[5] ]),
         var(weight[ feed==treat[6] ])) 
Var

# test of homogeneity of variances
# H0: sigma.1 = sigma.2 = sigma.3 = sigma.4 = sigma.5 = sigma.6 
# H1: there exist i,j s.t. sigma.i!=sigma.j
bartlett.test(weight, feed)


#### One-way ANOVA 
###-----------------
fit <- aov(weight ~ feed)
summary(fit)


#### How to read the summary:
#              Df   Sum Sq      Mean Sq      F value     Pr(>F)    
#  treat      (g-1) SStreat  SStreat/(g-1)  Fstatistic  p-value [H0: tau.i=0 for every i]
#  Residuals  (n-g) SSres     SSres/(n-g)                    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
###


#### 2 BY 2 COMPARISONS (g(g-1))/2 BONFERRONI
###-----------------
k <- g*(g-1)/2
alpha= 0.05
Mediag  <- tapply(d$weight, d$feed, mean)
SSres <- sum(residuals(fit)^2)
S <- SSres/(n-g)

# Example: CI for the difference "1 - 2"
paste(treat[1],"-",treat[2])
as.numeric(c(Mediag[1]-Mediag[2] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[1] + 1/ng[2] )),
             Mediag[1]-Mediag[2] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[1] + 1/ng[2] ))))


# CI for all the differences
ICrange=NULL
for(i in 1:(g-1)) {
  for(j in (i+1):g) {
    print(paste(treat[i],"-",treat[j]))        
    print(as.numeric(c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                       Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
    ICrange=rbind(ICrange,as.numeric(c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                                       Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
  }}


# Plot all differences
x11(width = 14, height = 7)
par(mfrow=c(1,2))
plot(feed, weight, xlab='treat', ylab='weight', col = rainbow(6), las=2)

h <- 1
plot(c(1,g*(g-1)/2),range(ICrange), pch='',xlab='pairs treat', ylab='Conf. Int. tau weight')
for(i in 1:(g-1)) {
  for(j in (i+1):g) {
    ind <- (i-1)*g-i*(i-1)/2+(j-i)
    lines (c(h,h), c(ICrange[ind,1],ICrange[ind,2]), col='grey55'); 
    points(h, Mediag[i]-Mediag[j], pch=16, col='grey55'); 
    points(h, ICrange[ind,1], col=rainbow(6)[j], pch=16);  # col = number of levels of g
    points(h, ICrange[ind,2], col=rainbow(6)[i], pch=16);  # Change col = rainbow(6)
    h <- h+1
  }}
abline(h=0)














################# One-way MANOVA ##############
# p=4, g=3

### Model: X.ij = mu + tau.i + eps.ij; eps.ij~N_p(0,Sigma), X.ij, mu, tau.i in R^p
### Test:
### H0: tau.1 = tau.2 = tau.3  = (0,0,0,0)'
### H1: (H0)^c
### that is
### H0: The membership to an iris species hasn't any significant effect on the mean
###     of X.ij (in any direction of R^4) 
### H1: There exists at least one direction in R^4 along which at least two species
###     have some feature significantly different

species.name <- factor(Species, labels=c('setosa','versicolor','virginica'))
iris4        <- iris[,1:4]

i1 <- which(species.name=='setosa')
i2 <- which(species.name=='versicolor')

x11(width=13)
par(mfrow=c(1,2))
boxplot(iris4[i1,], main='SETOSA',     ylim=c(0,8), col = rainbow(4))
boxplot(iris4[i2,], main='VERSICOLOR', ylim=c(0,8), col = rainbow(4))

n1 <- length(i1)
n2 <- length(i2)
n  <- n1+n2
g  <- length(levels(species.name))
p  <- 4



#### Verify assumptions:
###-----------------
# 1)  normality (multivariate) in each group (g tests)
Ps <- NULL
for(i in 1:g)
  Ps <- c(Ps, mcshapiro.test(iris[get(paste('i',i, sep='')),1:4])$p)  ## 1:4 = p number of numeric columns
Ps
# 2) same covariance structure (= same covariance matrix Sigma)
S  <-  cov(iris4)
S1 <-  cov(iris4[i1,])
S2 <-  cov(iris4[i2,])
# Qualitatively:
round(S1,digits=1)
round(S2,digits=1)
x11(width=21)
par(mfrow=c(1,3))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))



#### One-way MANOVA 
###----------------
fit <- manova(as.matrix(iris4) ~ species.name)
summary.manova(fit,test="Wilks") # Exact tests for p<=2 or g<=3 already implemented in R

### Who's the responsible for rejection this?
### Via ANOVA: for each of the p=4 variables we perform an ANOVA test
###            to verify if the membership to a group has influence
###            on the mean of the variable (we explore separately the
###            4 axes directions in R^4)
summary.aov(fit)


#### Bonferroni Confidence intervals on the TAU_i-TAU_j
# for each group for each of the p=4 numerical variables
###----------------
alpha <- 0.05
k <- p*g*(g-1)/2
qT <- qt(1-alpha/(2*k), n-g)

W <- summary.manova(fit)$SS$Residuals
m  <- sapply(iris4,mean)         # estimates mu
m1 <- sapply(iris4[i1,],mean)    # estimates mu.1=mu+tau.1
m2 <- sapply(iris4[i2,],mean)    # estimates mu.2=mu+tau.2
m3 <- sapply(iris4[i3,],mean)    # estimates mu.3=mu+tau.3

inf12 <- m1-m2 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
sup12 <- m1-m2 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
inf13 <- m1-m3 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n3) )
sup13 <- m1-m3 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n3) )
inf23 <- m2-m3 - qT * sqrt( diag(W)/(n-g) * (1/n2+1/n3) )
sup23 <- m2-m3 + qT * sqrt( diag(W)/(n-g) * (1/n2+1/n3) )

### CI <- list (tau1, tau2, tau3) each of 4 components
CI <- list(setosa_versicolor=cbind(inf12, sup12), setosa_virginica=cbind(inf13, sup13), versicolor_virginica=cbind(inf23, sup23))
CI


### Analog more generic
ng <- c(length(which(type==levels(factor(type))[1])),length(which(type==levels(factor(type))[2])),length(which(type==levels(factor(type))[3])))
n <- sum(ng)
alpha <- 0.1
k <- p*g*(g-1)/2
qT <- qt(1-alpha/(2*k), n-g)
W <- summary.manova(fit)$SS$Residuals
m  <- sapply(data,mean)         # estimates mu
m1 <- sapply(data[i1,],mean)    # estimates mu.1=mu+tau.1
m2 <- sapply(data[i2,],mean)    # estimates mu.2=mu+tau.2
inf12 <- m1-m2 - qT * sqrt( diag(W)/(n-g) * (1/ng[1]+1/ng[2]) )
sup12 <- m1-m2 + qT * sqrt( diag(W)/(n-g) * (1/ng[1]+1/ng[2]) )
CI <- list(cbind(inf12, sup12))
CI



















################# Two-ways ANOVA ##############
##### (p=1, g=2, b=2)
#####----------------------

km          <- c(18.7, 16.8, 20.1, 22.4, 14.0, 15.2, 22.0, 23.3)
distr       <- factor(c('Esso','Esso','Esso','Esso','Shell','Shell','Shell','Shell'))
benz        <- factor(c('95','95','98','98','95','95','98','98'))
distr_benz  <- factor(c('Esso95','Esso95','Esso98','Esso98','Shell95','Shell95','Shell98','Shell98'))


g <- length(levels(distr)) # factor 1
b <- length(levels(benz)) # factor 2
n <- length(km)/(g*b)

M           <- mean(km)
Mdistr      <- tapply(km,      distr, mean)
Mbenz       <- tapply(km,       benz, mean)
Mdistr_benz <- tapply(km, distr_benz, mean)


#### Two-ways ANOVA
###----------------
### Model with interaction (complete model): 
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N(0,sigma^2), 
###     i=1,2 (effect station), j=1,2 (effect gasoline)
fit.aov2.int <- aov(km ~ distr + benz + distr:benz)
summary.aov(fit.aov2.int)

### Test:
### 1) H0: gamma.11 = gamma.12 = gamma.21 = gamma.22 = 0    vs   H1: (H0)^c
###    i.e.,
###    H0: There is no significant interaction between the factors station
###        and gasoline in terms of performances
### 2) H0: tau.1 = tau.2 = 0    vs   H1: (H0)^c
###    i.e.,
###    H0: The effect "gas station" doesn't significantly influence performances 
### 3) H0: beta.1 = beta.2 = 0    vs   H1: (H0)^c
###    i.e.,
###    H0: The effect "gasoline" doesn't significantly influence performances

### Additive model: 
### X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N(0,sigma^2), 
###     i=1,2 (effect station), j=1,2 (effect gasoline)
fit.aov2.ad <- aov(km ~ distr + benz)
summary.aov(fit.aov2.ad)


#### Global test for significance of the two treatments
###----------------
SSdistr <- sum(n*b*(Mdistr - M)^2)              # or from the summary: 1.53    
SSbenz  <- sum(n*g*(Mbenz  - M)^2)              # or from the summary: 66.70
SSres   <- sum((km - M)^2) - (SSdistr+SSbenz)   # or from the summary: 16.37
Ftot      <- ( (SSdistr + SSbenz) / ((g-1)+(b-1)))/(SSres / (n*g*b-g-b+1))
Ptot      <- 1 - pf(Ftot, (g-1)+(b-1), n*g*b-g-b+1) # attention to the dof!
Ptot


### Interval at 90% for the differences (reduced additive model)
### [b=2, thus one interval only]
IC <- c(diff(Mbenz) - qt(0.95, (n*g-1)*b) * sqrt(SSres/((n*g-1)*b) *(1/(n*g) + 1/(n*g))), 
        diff(Mbenz) + qt(0.95, (n*g-1)*b) * sqrt(SSres/((n*g-1)*b) *(1/(n*g) + 1/(n*g))))
names(IC) <- c('Inf', 'Sup')
IC    # IC for mu(98)-mu(95)


##### THESE HYPOTHESES ARE REDUCED !!!NO!!!
### Note: we should have verified the hypotheses of normality and variance
###       homogeneity for the complete model, but with only 2 data for each
###       group we can't perform the tests.
### => we verify the assumptions on the reduced model (one-way ANOVA)
# 1) normality (univariate) in each group (2 tests)
Ps <- c(shapiro.test(km[ benz==levels(benz)[1] ])$p,
        shapiro.test(km[ benz==levels(benz)[2] ])$p)
Ps
# 2) homogeneity of variances
bartlett.test(km, benz)




















################# Two-ways MANOVA ##############
##### (p=3, g=2, b=2)
#####----------------------
### Model with interaction (complete model): 
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N_p(0,Sigma), [p=3]
###     i=1,2 (effect Extrusion), j=1,2 (effect Additive),
###     X.ijs, mu, tau.i, beta.j, gamma.ij in R^3

Ex   <- factor(plastic$Ex, labels=c('L','H')) # Treat.1
Ad   <- factor(plastic$Ad, labels=c('L','H')) # Treat.2
ExAd <- Ex
levels(ExAd) <- c('LL','LH','HL','HH')
ExAd[Ex=='L' & Ad=='L'] <- 'LL'
ExAd[Ex=='L' & Ad=='H'] <- 'LH'
ExAd[Ex=='H' & Ad=='L'] <- 'HL'
ExAd[Ex=='H' & Ad=='H'] <- 'HH'
plastic3  <- plastic[,3:5]

#### Verify the assumptions (although we only have 5 data in each group!)
####----------------------
# 1) normality (multivariate) in each group (4 test)
Ps <- c(mcshapiro.test(plastic3[ ExAd==levels(ExAd)[1], ])$p,
        mcshapiro.test(plastic3[ ExAd==levels(ExAd)[2], ])$p,
        mcshapiro.test(plastic3[ ExAd==levels(ExAd)[3], ])$p,
        mcshapiro.test(plastic3[ ExAd==levels(ExAd)[4], ])$p)
Ps

# 2) homogeneity of the covariance (qualitatively)
S1 <-  cov(plastic3[ ExAd==levels(ExAd)[1], ])
S2 <-  cov(plastic3[ ExAd==levels(ExAd)[2], ])
S3 <-  cov(plastic3[ ExAd==levels(ExAd)[3], ])
S4 <-  cov(plastic3[ ExAd==levels(ExAd)[4], ])
x11(width=21)
par(mfrow=c(1,4))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3,S4), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3,S4), (0:100)/100, na.rm=TRUE))
image(S3, col=heat.colors(100),main='Cov. S3', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3,S4), (0:100)/100, na.rm=TRUE))
image(S4, col=heat.colors(100),main='Cov. S4', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3,S4), (0:100)/100, na.rm=TRUE))
bartlett.test(plastic3,ExAD)



#### Two-ways MANOVA 
###-----------------
### Model with interaction (complete model): 
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N_p(0,Sigma), [p=3]
###     i=1,2 (effect Extrusion), j=1,2 (effect Additive),
###     X.ijs, mu, tau.i, beta.j, gamma.ij in R^3
fit <- manova( as.matrix(plastic3) ~ Ex + Ad + Ex:Ad)
summary.manova(fit, test="Wilks")

### Model without interaction (additive model): 
### X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N_p(0,Sigma), [p=3]
###     i=1,2 (effect Extrusion), j=1,2 (effect additive),
###     X.ijs, mu, tau.i, beta.j, in R^3
fit2<- manova( as.matrix(plastic3) ~ Ex + Ad)
summary.manova(fit2, test="Wilks")
summary.aov(fit2)





#### BONFERRONI CONF. FOR DIFFERENCE OF MEANS
###-----------------
alpha <- 0.05
g <- 2
b <- 2
p <- 3
n <- 5
N <- n*g*b # 20
W <- summary.manova(fit2)$SS$Residuals

# how many comparisons?
k <- g*(g-1)/2*p + b*(b-1)/2*p
# because we have: g levels on the first treatment on p components
#                  b levels on the second treatment on p components
k
qT <- qt(1 - alpha / (2 * k), g*b*n-g-b+1)
# the degrees of freedom of the residuals on the additive model are
# g*b*n-g-b+1
mExL  <- sapply(plastic3[Ex=='L',],mean)
mExH  <- sapply(plastic3[Ex=='H',],mean)
infEx <- mExH-mExL - qT * sqrt( diag(W)/(g*b*n-g-b+1) * (1/10+1/10) )
supEx <- mExH-mExL + qT * sqrt( diag(W)/(g*b*n-g-b+1) * (1/10+1/10) )

mAdL  <- sapply(plastic3[Ad=='L',],mean)
mAdH  <- sapply(plastic3[Ad=='H',],mean)
infAd <- mAdH-mAdL - qT * sqrt( diag(W)/(g*b*n-g-b+1) * (1/10+1/10) )
supAd <- mAdH-mAdL + qT * sqrt( diag(W)/(g*b*n-g-b+1) * (1/10+1/10) )

IC2   <- list(ExH_ExL=cbind(infEx, supEx), AdH_AdL=cbind(infAd, supAd))
IC2








######## BONFERRONI CONF. INTERVALS FOR MEANS AND VARIANCE
fit3 <- aov(weight ~ species, penguins) 
summary(fit3)

DF <- fit3$df # n*g*b - 1 - (g-1) = 15*3*2-1-2 = 90-3 = 87
Spooled <- sum(fit3$res^2)/DF

means <- as.vector(tapply(penguins$weight, penguins$species, mean))
names(means) <- levels(species)
means

alpha <- 0.10
k     <- 4 # g + 1 = 4 (g Conf Int for the means and 1 for the variance)
BF    <- rbind(cbind(means - sqrt(Spooled / 30) * qt(1 - alpha / (2*k), DF), 
                     means + sqrt(Spooled / 30) * qt(1 - alpha / (2*k), DF)),
               c(Spooled * DF / qchisq(1 - alpha / (2*k), DF), 
                 Spooled * DF / qchisq(alpha / (2*k), DF)))
rownames(BF)[4] <- 'Var.'
BF