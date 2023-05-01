#########
## One-way ANOVA - one variable (weight) observed over g levels (feed)
#########

chickwts <- read.table(file='chickwts.txt', header=T)
head(chickwts)
dim(chickwts)

# BOXPLOT
plot(feed, weight, xlab='treat', ylab='weight', col='grey85', main='Dataset Chicken Weights')


### Model: weigth.ij = mu + tau.i + eps.ij; eps.ij~N(0,sigma^2)
### Test:
### H0: tau.1 = tau.2 = tau.3 = tau.4 = tau.5 = tau.6 = 0
### H1: (H0)^c
### i.e.,
### H0: The feed supplements don't have effect (= "chickens belong to a single population")
### H1: At least one feed supplement has effect (= "chickens belong to 2, 3, 4, 5 or 6 populations")

# PLOT OF THE TWO MODELS
par(mfrow=c(1,2))
barplot(rep(mean(weight),6), names.arg=levels(feed), ylim=c(0,max(weight)),
        las=2, col='grey85', main='Model under H0')
barplot(tapply(weight, feed, mean), names.arg=levels(feed), ylim=c(0,max(weight)),
        las=2, col=rainbow(6),main='Model under H1')

# DATA
n       <- length(feed)      # total number of obs.
ng      <- table(feed)       # number of obs. in each group
treat   <- levels(feed)      # levels of the treatment
g       <- length(treat)     # number of levels (i.e., of groups)

### verify the assumptions:
# 1) normality (univariate) in each group (6 tests)
Ps <- c(shapiro.test(weight[ feed==treat[1] ])$p,
        shapiro.test(weight[ feed==treat[2] ])$p,
        shapiro.test(weight[ feed==treat[3] ])$p,
        shapiro.test(weight[ feed==treat[4] ])$p,
        shapiro.test(weight[ feed==treat[5] ])$p,
        shapiro.test(weight[ feed==treat[6] ])$p) 
Ps

# 2) same covariance structure (= same sigma^2)
# test of homogeneity of variances
# H0: sigma.1 = sigma.2 = sigma.3 = sigma.4 = sigma.5 = sigma.6 
# H1: there exist i,j s.t. sigma.i!=sigma.j
bartlett.test(weight, feed)    #bartlett.test(dou$prezzo, factor(dou$citta):factor(dou$tipo))

#One_way ANOVA
fit <- aov(weight ~ feed)
summary(fit)

### How to read the summary:
#              Df   Sum Sq      Mean Sq      F value     Pr(>F)    
#  treat      (g-1) SStreat  SStreat/(g-1)  Fstatistic  p-value [H0: tau.i=0 for every i]
#  Residuals  (n-g) SSres     SSres/(n-g)                    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
###

### Which supplement is responsible for this? Use Bonferroni
k <- g*(g-1)/2
alpha= 0.05

Mediag  <- tapply(weight, feed, mean)
SSres <- sum(residuals(fit)^2)
S <- SSres/(n-g)


# CI for the difference "casein [1] - horsebean [2]"


paste(treat[1],"-",treat[2])
as.numeric(c(Mediag[1]-Mediag[2] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[1] + 1/ng[2] )),
             Mediag[1]-Mediag[2] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[1] + 1/ng[2] ))))




# CI for all the differences

ICrange = NULL
for(i in 1:(g-1)) {
    for(j in (i+1):g) {
        print(paste(treat[i],"-",treat[j]))
        interval = c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                     Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )))
        interval = as.numeric(interval)
        print(interval)
        ICrange = rbind(ICrange,interval)
    }
}

# ---------------------------------------------------------------------------



#########
## Two-ways ANOVA - one variable (distance) observed over g * b levels (gas station & gasoline)
#########

km          <- c(18.7, 16.8, 20.1, 22.4, 14.0, 15.2, 22.0, 23.3)
distr       <- factor(c('Esso','Esso','Esso','Esso','Shell','Shell','Shell','Shell'))
benz        <- factor(c('95','95','98','98','95','95','98','98'))
distr_benz  <- factor(c('Esso95','Esso95','Esso98','Esso98','Shell95','Shell95','Shell98','Shell98'))

g <- length(levels(distr))
b <- length(levels(benz))
n <- length(km)/(g*b)

M           <- mean(km)
Mdistr      <- tapply(km,      distr, mean)
Mbenz       <- tapply(km,       benz, mean)
Mdistr_benz <- tapply(km, distr_benz, mean)

############## if many data, verify assumptions -> g*b shapiro, 1 bartlett

### Model with interaction (complete model): 
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N(0,sigma^2), 
###     i=1,2 (effect station), j=1,2 (effect gasoline)

fit.aov2.int <- aov(km ~ distr + benz + distr:benz)
summary.aov(fit.aov2.int)

### Test:
### 1) H0: gamma.11 = gamma.12 = gamma.21 = gamma.22 = 0    vs   H1: (H0)^c
# Focus on distr:benz, pvalue=0.05857  -> Reject at 10%, don't reject at 1%,5%

### 2) H0: tau.1 = tau.2 = 0    vs   H1: (H0)^c
# Focus on distr, pvalue=0.37001. Don't reject at 10%, 5%, 1% -> not significant

### 3) H0: beta.1 = beta.2 = 0    vs   H1: (H0)^c
# Focus on benz, pvalue=0.00264. Reject at 10%, 5%, 1% -> significant

#remove the interaction 

## Additive model: 
### X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N(0,sigma^2), 
###     i=1,2 (effect station), j=1,2 (effect gasoline)
fit.aov2.ad <- aov(km ~ distr + benz)
summary.aov(fit.aov2.ad)

# Test: 2bis) H0: tau.1 = tau.2 = 0    vs   H1: (H0)^c
# Focus on distr, pvalue=0.52440. Don't reject at 10%, 5%, 1% -> not significant
# Test: 3bis) H0: beta.1 = beta.2 = 0    vs   H1: (H0)^c
# Focus on benz, pvalue0.00632. Reject at 10%, 5%, 1% -> significant

### Reduced additive model (ANOVA one-way, b=2): 
### X.jk = mu + beta.j + eps.jk; eps.jk~N(0,sigma^2), 
###     j=1,2 (effect gasoline)
fit.aov1 <- aov(km ~ benz)
summary.aov(fit.aov1)

SSres <- sum(residuals(fit.aov1)^2)

### Interval at 90% for the differences (reduced additive model) [b=2, thus one interval only]
IC <- c(diff(Mbenz) - qt(0.95, (n*g-1)*b) * sqrt(SSres/((n*g-1)*b) *(1/(n*g) + 1/(n*g))), 
        diff(Mbenz) + qt(0.95, (n*g-1)*b) * sqrt(SSres/((n*g-1)*b) *(1/(n*g) + 1/(n*g))))
names(IC) <- c('Inf', 'Sup')
IC    # IC for mu(98)-mu(95)

### => we verify the assumptions on the reduced model (one-way ANOVA)
# 1) normality (univariate) in each group (2 tests)
Ps <- c(shapiro.test(km[which(benz==levels(benz)[1])] )$p,
        shapiro.test(km[which(benz==levels(benz)[2])] )$p)
Ps

# 2) homogeneity of variances
bartlett.test(km, benz)    #bartlett.test(dou$prezzo, factor(dou$citta):factor(dou$tipo))

graphics.off()



#???????
#--------------------------------------------------

#CI on means and variance
g <- length(levels(factor(tipo)))
DF <- fit3$df # n*g*b - 1 - (g-1)  #n-1
Spooled <- sum(fit3$res^2)/DF      #(n-1)*diag(S)

means <- as.vector(tapply(penguins$weight, penguins$species, mean))
names(means) <- levels(species)

alpha <- 0.10
k     <- g + 1 

ng <- c(length(which(tipo==levels(factor(tipo))[1])),
        length(which(tipo==levels(factor(tipo))[2])),
        length(which(tipo==levels(factor(tipo))[3])))

qT <- qt(1-alpha/(2*k), DF)

BF    <- rbind(cbind(inf=means - sqrt(c(Spooled) / ng) * qT,
                     sup=means + sqrt(c(Spooled) / ng) * qT),
               c(Spooled * DF / qchisq(1 - alpha / (2*k), DF), 
                 Spooled * DF / qchisq(alpha / (2*k), DF)))


## PARAMETERS
m  <- sapply(iris4,mean)         # estimates mu
m1 <- sapply(iris4[i1,],mean)    # estimates mu.1=mu+tau.1
m2 <- sapply(iris4[i2,],mean)    # estimates mu.2=mu+tau.2
m3 <- sapply(iris4[i3,],mean)    # estimates mu.3=mu+tau.3

