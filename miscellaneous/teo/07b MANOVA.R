#########
## One-way MANOVA - p variables (lenght, width, ...) observed over g levels (species)
#########

iris <- read.table(file='iris.txt', header=T)
head(iris)
dim(iris)

species.name <- factor(Species, labels=c('setosa','versicolor','virginica'))
iris4        <- iris[,1:4]

i1 <- which(species.name=='setosa')
i2 <- which(species.name=='versicolor')
i3 <- which(species.name=='virginica')

### Model: X.ij = mu + tau.i + eps.ij; eps.ij~N_p(0,Sigma), X.ij, mu, tau.i in R^4
### Test:
### H0: tau.1 = tau.2 = tau.3  = (0,0,0)'
### H1: (H0)^c
### that is
### H0: The membership to an iris species hasn't any significant effect on the mean
###     of X.ij (in any direction of R^4) 
### H1: There exists at least one direction in R^4 along which at least two species
###     have some feature significantly different

n1 <- length(i1)
n2 <- length(i2)
n3 <- length(i3)
n  <- n1+n2+n3

g  <- length(levels(species.name))
p  <-  length(iris4)

### Verify the assumptions:
# 1)  normality (multivariate) in each group (3 tests)
Ps <- NULL
for(i in 1:g)
  Ps <- c(Ps, mcshapiro.test(iris[get(paste('i',i, sep='')),1:4])$p) 
Ps

# 2) same covariance structure (= same covariance matrix Sigma)
S  <-  cov(iris4)
S1 <-  cov(iris4[i1,])
S2 <-  cov(iris4[i2,])
S3 <-  cov(iris4[i3,])

# Qualitatively:

round(S1,digits=1)
round(S2,digits=1)
round(S3,digits=1)

x11(width=21)
par(mfrow=c(1,3))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S3, col=heat.colors(100),main='Cov. S3', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))

dev.off()

### One-way MANOVA 
fit <- manova(as.matrix(iris4) ~ species.name)
summary.manova(fit,test="Wilks")

### Who's the responsible for this?
# analyse each direction with an anova
summary.aov(fit)

### Via Bonferroni
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

CI <- list(setosa_versicolor=cbind(inf12, sup12), setosa_virginica=cbind(inf13, sup13), versicolor_virginica=cbind(inf23, sup23))
CI

#########
## Two-ways MANOVA - p variables (lenght, width, ...) observed over g * b levels (species) 
#########

#line 576 lab 7, here

#_______________________________________________________________________________
##### Two-ways MANOVA
##### (p=3, g=2, b=2)
#####----------------------

### Example 6.13 JW (p.319)
###------------------------
# The optimum conditions for extruding plastic films have been 
# examined using a technique called Evolutionary Operation. In 
# the course of the study, three responses were measured:
# X.1 = Tear Resistance
# X.2 = Gloss
# X.3 = Opacity
# at two levels of the factors:
# factor.1 = Rate of Extrusion     (0=Low level, 1=High level)
# factor.2 = Amount of an Additive (0=Low level, 1=High level)
# The measurements were repeated n=5 times at each combination 
# of the factor levels.

plastic <- read.table('T6-4.dat',col.names=c('Ex','Ad','Tr','Gl','Op'))
plastic <- read.table(file.choose(),col.names=c('Ex','Ad','Tr','Gl','Op'))
plastic

Ex   <- factor(plastic$Ex, labels=c('L','H')) # Treat.1
Ad   <- factor(plastic$Ad, labels=c('L','H')) # Treat.2

ExAd <- Ex
levels(ExAd) <- c('LL','LH','HL','HH')
ExAd[Ex=='L' & Ad=='L'] <- 'LL'
ExAd[Ex=='L' & Ad=='H'] <- 'LH'
ExAd[Ex=='H' & Ad=='L'] <- 'HL'
ExAd[Ex=='H' & Ad=='H'] <- 'HH'

plastic3  <- plastic[,3:5]

### Graphical exploration of the data
# effect of the treatments + their interaction on the first variable
x11()
layout(matrix(c(1,1,2,3), 2, byrow=T))
boxplot(plastic3[,1]~ExAd, main='Model with Interac. Extrusion+Additive (Tear Resistance)', ylab='Tr', col='grey95')
boxplot(plastic3[,1]~Ex,   main='Only Factor Extrusion'  , ylab='Tr', col=c('red','blue'))
boxplot(plastic3[,1]~Ad,   main='Only Factor Additive'   , ylab='Tr', col=c('forestgreen','gold'))

# effect of the treatments + their interaction on the second variable
x11()
layout(matrix(c(1,1,2,3), 2, byrow=T))
boxplot(plastic3[,2]~ExAd, main='Model with Interac. Extrusion+Additive (Gloss)', ylab='Gl', col='grey95')
boxplot(plastic3[,2]~Ex,   main='Only Factor Extrusion'  , ylab='Gl', col=c('red','blue'))
boxplot(plastic3[,2]~Ad,   main='Only Factor Additive'   , ylab='Gl', col=c('forestgreen','gold'))

# effect of the treatments + their interaction on the third variable
x11()
layout(matrix(c(1,1,2,3), 2, byrow=T))
boxplot(plastic3[,3]~ExAd, main='Model with Interac. Extrusion+Additive (Opacity)', ylab='Op', col='grey95')
boxplot(plastic3[,3]~Ex,   main='Only Factor Extrusion'  , ylab='Op', col=c('red','blue'))
boxplot(plastic3[,3]~Ad,   main='Only Factor Additive'   , ylab='Op', col=c('forestgreen','gold'))

dev.off()
dev.off()
dev.off()

### Model with interaction (complete model): 
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N_p(0,Sigma), [p=3]
###     i=1,2 (effect Extrusion), j=1,2 (effect Additive),
###     X.ijs, mu, tau.i, beta.j, gamma.ij in R^3

### Verify the assumptions (although we only have 5 data in each group!)
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

dev.off()

### Two-ways MANOVA 
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

# Both the treatments have a significant effect on the mean (but not
# their interaction, that we could remove)

# Let's verify if this is true for all the variables through appropriate
# conf. int. and tests on the components:

# ANOVA on the components (we look at the 3 axes-directions in R^3
#                          separately)
summary.aov(fit2)

# Bonferroni
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


x11(width=21, height = 14)
par(mfrow=c(3,4))
boxplot(plastic3[,1]~Ex,   main='Fact.: Extrusion (Tear Resistance)'  , ylab='Tr', col=rainbow(2*6)[c(1,2)], ylim=c(-2,10))
plot(c(1,g*(g-1)/2),range(IC2[[1]][1,]), pch='',main='IC (tau.1-tau.2)[1]',xlab='pairs treat', ylab='IC (tau.1-tau.2)[1]', ylim=c(-2,10))
lines(c(1,1), c(IC2[[1]][1,1],IC2[[1]][1,2]), col='grey55'); 
points(1, (mExH-mExL)[1], pch=16, col='grey55'); 
points(1, IC2[[1]][1,1], col=rainbow(2*6)[1], pch=16); 
points(1, IC2[[1]][1,2], col=rainbow(2*6)[2], pch=16);
abline(h=0)

boxplot(plastic3[,1]~Ad,   main='Fact.: Additive (Tear Resistance)'  , ylab='Tr', col=rainbow(2*6)[c(7,8)], ylim=c(-2,10))
plot(c(1,g*(g-1)/2),range(IC2[[2]][1,]), pch='',main='IC (beta.1-beta.2)[1]',xlab='pairs treat', ylab='IC (beta.1-beta.2)[1]', ylim=c(-2,10))
lines(c(1,1), c(IC2[[2]][1,1],IC2[[2]][1,2]), col='grey55'); 
points(1, (mAdH-mAdL)[1], pch=16, col='grey55'); 
points(1, IC2[[2]][1,1], col=rainbow(2*6)[7], pch=16); 
points(1, IC2[[2]][1,2], col=rainbow(2*6)[8], pch=16);
abline(h=0)

boxplot(plastic3[,2]~Ex,   main='Fact.: Extrusion (Gloss)'  , ylab='Gl', col=rainbow(2*6)[c(3,4)], ylim=c(-2,10))
plot(c(1,g*(g-1)/2),range(IC2[[1]][2,]), pch='',main='IC (tau.1-tau.2)[2]',xlab='pairs treat', ylab='IC (tau.1-tau.2)[2]', ylim=c(-2,10))
lines(c(1,1), c(IC2[[1]][2,1],IC2[[1]][2,2]), col='grey55'); 
points(1, (mExH-mExL)[2], pch=16, col='grey55'); 
points(1, IC2[[1]][2,1], col=rainbow(2*6)[3], pch=16); 
points(1, IC2[[1]][2,2], col=rainbow(2*6)[4], pch=16);
abline(h=0)

boxplot(plastic3[,2]~Ex,   main='Fact.: Additive (Gloss)'  , ylab='Gl', col=rainbow(2*6)[c(9,10)], ylim=c(-2,10))
plot(c(1,g*(g-1)/2),range(IC2[[2]][2,]), pch='',main='IC (beta.1-beta.2)[2]',xlab='pairs treat', ylab='IC (beta.1-beta.2)[2]', ylim=c(-2,10))
lines(c(1,1), c(IC2[[2]][2,1],IC2[[2]][2,2]), col='grey55'); 
points(1, (mAdH-mAdL)[2], pch=16, col='grey55'); 
points(1, IC2[[2]][2,1], col=rainbow(2*6)[9], pch=16); 
points(1, IC2[[2]][2,2], col=rainbow(2*6)[10], pch=16);
abline(h=0)

boxplot(plastic3[,3]~Ex,   main='Fact.: Extrusion (Opacity)'  , ylab='Op', col=rainbow(2*6)[c(5,6)], ylim=c(-2,10))
plot(c(1,g*(g-1)/2),range(IC2[[1]][3,]), pch='',main='IC (tau.1-tau.2)[3]',xlab='pairs treat', ylab='IC (tau.1-tau.2)[3]', ylim=c(-2,10))
lines(c(1,1), c(IC2[[1]][3,1],IC2[[1]][3,2]), col='grey55'); 
points(1, (mExH-mExL)[3], pch=16, col='grey55'); 
points(1, IC2[[1]][3,1], col=rainbow(2*6)[5], pch=16); 
points(1, IC2[[1]][3,2], col=rainbow(2*6)[6], pch=16);
abline(h=0)

boxplot(plastic3[,3]~Ex,   main='Only Factor Additive (Opacity)'  , ylab='Op', col=rainbow(2*6)[c(11,12)], ylim=c(-2,10))
plot(c(1,g*(g-1)/2),range(IC2[[2]][3,]), pch='',main='IC (beta.1-beta.2)[3]',xlab='pairs treat', ylab='IC beta.1[3]', ylim=c(-2,10))
lines(c(1,1), c(IC2[[2]][3,1],IC2[[2]][3,2]), col='grey55'); 
points(1, (mAdH-mAdL)[3], pch=16, col='grey55'); 
points(1, IC2[[2]][3,1], col=rainbow(2*6)[11], pch=16); 
points(1, IC2[[2]][3,2], col=rainbow(2*6)[12], pch=16);
abline(h=0)

dev.off()





