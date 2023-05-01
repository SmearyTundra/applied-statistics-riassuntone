load("/Users/Nicole/OneDrive - Politecnico di Milano/Corsi/AS/Laboratori/mcshapiro.test.RData")
data <- read.table("data.txt",header = T)
par(mfrow=c(1,3))

# SALVARE GRAFICI 
pdf(file = "myplot2.pdf", onefile = T)
plot(rnorm(10),rnorm(10))
plot(rnorm(10),rnorm(10))
plot(rnorm(10),rnorm(10))
dev.off()

###LAB-1: DATA IMPORT/EXPORT ----
record <- read.table('record.txt', header=T)
record

head(record)
dim(record)
dimnames(record)

write.table(record, file = 'record_mod.txt') #save data frame

#Vettori
u <- seq(2,5,by=1)#sequenza da 2 a 5,con incremento di 1 or z <- 2:5

## FACTOR
district <- c('MI',  'MI',  'VA',  'BG',  'LO', 'LO',  'CR',  'Alt', 'CR',  'MI',  
              'Alt', 'CR',  'LO',  'VA',  'MI',  'Alt', 'LO',  'MI')
district <- factor(district,levels=c('MI','LO','BG','CR','VA','Alt'))
district

## ANALISI QUANTITATIVA 
colMeans(record)#media per colonne
sapply(record, mean)#media per colonne
sapply(record, sd)
sapply(record, var)
cov(record)#matrice cov
cor(record)#matrice cor
 
#Grafici
x11()
par(mfrow=c(2,2))
hist(covariata1,prob=T,main="Histogram records 100m",xlab="sec")
hist(covariata2,prob=T,main="Histogram records 200m",xlab="sec")
boxplot(record[,1:2],main="Boxplot records covariata1 e covariata2",xlab="sec")
plot(m100,m200, main='Scatter plot records covariata1 e covariata2',xlab="Records covariata1",ylab="Records covariata2")

#dati gaussiani
qqnorm(covariata) # quantile-quantile plot
qqline(covariata, col='red')
shapiro.test(covariata1)

##### Visualization of Multivariate Data
#Scatterplot
pairs(record)
#Boxplot
boxplot(record,col='gold')
boxplot(log(record), col='gold')
# Starplot
stars(record, col.stars=rep('gold',55))
# Radarplot
stars(record, draw.segments=T)

#save plot
pdf(file = "myplot2.pdf", onefile = T)
plot(rnorm(10),rnorm(10))
plot(rnorm(10),rnorm(10))
plot(rnorm(10),rnorm(10))

#UNIVARIATE TEST
T.x <- abs (mean(x))/sqrt(var(x/n))#test
P <- (1-pt(T.x,n-1))*2#è-value
#or
t.test(x, alternative="two.side", mu = 0, conf.level = 0.9)

color.position <- ifelse(aneurysm.position == '1', 'red', 'blue')


x11()
plot(1:8,sort(p_val))# metto in rodine crescente i p_val 
x=seq(1,8,by=0.25)#creo retta con pendenza data da y 
y=x*0.25/8#alpha/k numero di p_value
lines(x,y)
#not reject for 6,7,8
#considero statisticamente migliorati quelli che sono sotto la linea

###LAB-3: PCA ----
# Scatter plot 
pairs(dataset, col=rainbow(dim(dataset)[1]), pch=16, main='Scatter plot')
#Boxplot
boxplot(dataset, las=2, col='gold')
boxplot(scale(x=dataset,center = T, scale=F), las=2, col='gold') #media=0

S <- cov(dataset)
image(S, asp=1)#matrice correlazione
var.gen <- det(S)
var.tot <- sum( diag(S) )

#pca --> trovo componenti utili
pc.dataset <- princomp(dataset, scores=T)
summary(pc.dataset)
#grafico loading
load.tour <- pc.tourists$loadings
load.tour[,1:x]#x--> n* di comp.che tengo

x11()
par(mfcol = c(4,2))
for(i in 1:x) barplot(load.tour[,i], ylim = c(-1, 1), main=paste("PC",i))

#Grafico varianza
x11()
layout(matrix(c(2,3,1,3),2,byrow=T))
plot(pc.dataset, las=2, main='Principal components', ylim=c(0,4.5e7))
barplot(sapply(dataset,sd)^2, las=2, main='Original Variables', ylim=c(0,4.5e7), ylab='Variances')
plot(cumsum(pc.dataset$sd^2)/sum(pc.dataset$sd^2), type='b', axes=F, xlab='number of components', 
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(dataset),labels=1:ncol(dataset),las=2)

#Standardize
tourists.sd <- scale(tourists)
tourists.sd <- data.frame(tourists.sd)

#Grafico per vedere dove si sviluppano le direzioni
biplot(pc.tourists)

# scores
scores.tourists <- pc.tourists$scores
scores.tourists

x11()
plot(scores.tourists[,1:2])
text(scores.tourists[,1:2],dimnames(dataset)[[1]], cex=0.7)
abline(h=0, v=0, lty=2, col='grey')

#scores on new data data.3jul
scores.3jul <- t(pca.NO$loadings)%*%(c-colMeans(NO))

#Matplot: uso quando c'è sequenza temporale
matplot(t(dataset), type='l', xlab='Age', ylab='Number of Residents', lty=1, col=rainbow(33), las=1)

# Project the new datum on the reduced space identified at point (c).

z0 <- c(MeanTemp=30,MinTemp=23,MaxTemp=36,DewPoint=22,Humidity=65, Visibility=19, MeanWind=5, MaxWind=15)
z0 <-scale(z0)
new_sc <- z0[,1] %*% load.dataset[,1:2]
new_sc

x11()
plot(scores.weather[,1:2])
points(new_sc[1],new_sc[2],col='red')
abline(h=0, v=0, lty=2, col='grey')

###LAB-4: Testing for multivariate normality ----
library(mvtnorm)
library(mvnormtest)

#Controllo graficamente normalità
hist(cov1, prob=T, ylab='density', xlab='X.1', main='Histogram of X.1',ylim=c(0,0.45))
lines((-1000):1000 /100, dnorm((-1000):1000 /100,mean(X[,1]),sd(X[,1])), col='blue', lty=2)

qqnorm(cov1, main='QQplot of X.1',xlab='theoretical quantiles', ylab='sample quantiles')
qqline(cov1)
#--> posso rifare il procediemento con le PC
shapiro.test(cov1)
mcshapiro.test(cov1)

# Data non gaussiani:
#1) Identify clusters 
#2)Identify (and possibly remove) outliers
#3)Transform the data (e.g., Box-Cox transformations, see Johnson-Wichern Chap.4.8, R functions powerTransform(); bcPower())
#4)Work without the Gaussian assumption (e.g., permutation tests)


###LAB-5: Box-Cox transformations ----
library(car)
# x_lambda = (x^lambda-1)/lambda if lambda!=0
#            ln(x)               if lambda==0
# For lambda<1: observations <1 are "spread", observations >1 are "shrinked"
# For lambda>1: observations <1 are "shrinked", observations >1 are "spread"
lambda.x <- powerTransform(dataset) #optimal lambda
bc.x <- bcPower(dataset, lambda.x$lambda) #transform sample with optimal lambda
#Bivariate Box-cox
lambda <- powerTransform(cbind(x,y))
BC.x <- bcPower(x, lambda$lambda[1])
BC.y <- bcPower(y, lambda$lambda[2])
#se lambda ~ 1 : lascio uguale, lambda ~ 0: log(cov)


###LAB-5: Tests and confidence regions for the mean of a multivariate Gaussian ---- 
mcshapiro.test(x) #verifico hp gauss

n <- dim(x)[1]
p <- dim(x)[2]
alpha <- 0.01
mu0 <- c(1,0)
x.mean   <- colMeans(x)
x.cov    <- cov(x)
x.invcov <- solve(x.cov)
#Verify the gaussianity assumption (mcshapiro)

### Test on the mean of level alpha=1%
### H0: mu == mu0 vs H1: mu != mu0,with mu0=c(0,0)
x.T2 <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) #T2 Statistics--> Hotteling's thm
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p) #Radius of the ellipsoid
# Rejection region: {x.T2>cfr.fisher}:
x.T2 < cfr.fisher # if TRUE : non rigetto H0
P <- 1-pf(x.T2*(n-p)/((n-1)*p), p, n-p) #rifiuto al livello alpha (> p.value)
#esemopio: pvalue=0.0234 --> rigetto al livello 5%,3% ma NON 1% e 2%

# Region of rejection (centered in mu0)
x11()
plot(x, asp = 1)
ellipse(mu0, shape=x.cov/n, sqrt(cfr.fisher), col = 'blue', lty = 2, center.pch = 16)

points(x.mean[1], x.mean[2], pch = 16, col ='red', cex = 1.5)# We add a red point in correspondence of the sample mean

# Confidence region (centered in x.mean)
# { m \in R^2 s.t. n * (x.mean-m)' %*% (x.cov)^-1 %*% (x.mean-m) < cfr.fisher }
ellipse(x.mean, x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 2, lwd=2, center.cex=1)
#plot3d(ellipse3d(centre=x.mean, x=x.cov/n, t=sqrt(cfr.fisher),level=0.95),col="blue",aspect=TRUE)

# Remark: the radius and the shape of the ellipse are the same, but the center changes:
# - Rejection region: the center is the mean mu0 under H0 (blue ellipse) (outside the ellipse)
# - Confidence region: the center is the sample mean (red ellipse)

# Which relation between the two ellipses?
# - If the rejection region does NOT contain the sample mean (i.e., we
#   are in the acceptance region), then we cannot reject H0 
#   (i.e., if the sample mean falls within the ellipse we accept H0)
# - If the mean under H0 (mu0) is contained in the confidence region
#   of level 1-alpha, then we do not reject H0 at level alpha
# => the confidence region of level 1-alpha contains all the mu0
#    that we would accept at level alpha

# Center:
x.mean

# Directions of the principal axes:
eigen(x.cov/n)$vectors
#Radius
r <- sqrt(cfr.fisher)
# Length of the semi-axes of the ellipse:
r*sqrt(eigen(x.cov/n)$values) 

#Asymptotic test (non verifico gauss)
### H0: mu == mu0 vs H1: mu != mu0
### with mu0=c(1,0)
x.T2A   <- n * (x.mean-mu0) %*%  x.invcov  %*% (x.mean-mu0)
cfr.chisq <- qchisq(1-alpha,p)
x.T2A < cfr.chisq 
PA <- 1-pchisq(x.T2A, p) #pvalue

#Simultaneous T2 confidence intervals on the coordinate directions
T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
# Both the intervals contain the mean under H0 --> non reject but ONE_AT_A_TIME
# (i.e., mu0 is contained in the rectangular region determined by
# the projection of the ellipsoid along the coordinate directions)

# Remark: this is not in contrast with the previous findings
# Rejecting the global T2-test means that we reject H0 along at least one
# direction, not necessarily along the coordinate direction

#intervalli grafici (anche per dim>2)
x11()
matplot(1:4,1:4,pch='',ylim=range(database),xlab='Variables',ylab='T2 for a component', 
        main='Simultaneous T2 conf. int. for the components')
for(i in 1:4)segments(i,T2[i,1],i,T2[i,3],lwd=3,col=i)
points(1:4, T2[,2], pch=16, col=1:4)

#Bonferroni intervals
k <- p # number of intervals I want to compute (set in advance)
cfr.t <- qt(1-alpha/(2*k),n-1)
Bf.mean <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
# Both the intervals contain the mean under H0
# (i.e., mu0 is contained in the rectangular region determined by
# the Bonferroni intervals along the coordinate directions)
# Remark: if we wanted to compute additional Bonferroni intervals
# along other directions, we would need to re-compute all the Bonferroni
# intervals with another correction k

#Bonferroni for the variance
Bf.var <- cbind(inf=diag(x.cov)*(n-1) / qchisq(1 - alpha/(2*k), n-1),
               center=diag(x.cov),
               sup=diag(x.cov)*(n-1) / qchisq(alpha/(2*k), n-1))

C <- rbind(c(1,0),c(0,1),c(1,1))
T2 <- cbind( 
  C%*%x.mean - sqrt(diag(C%*%x.cov%*%t(C))/n*(p*(n-1)/(n-p))*qf(1-alpha,p,n-p)),
  C%*%x.mean ,
  C%*%x.mean + sqrt(diag(C%*%x.cov%*%t(C))/n*(p*(n-1)/(n-p))*qf(1-alpha,p,n-p)))
colnames(T2) <- c('Inf','Mean','Sup')
T2



###LAB-6: Test for the mean of paired multivariate Gaussian observations ----
library(car)
# compute the sample of differences of measurement

# Test the Gaussian assumption (on D!)
mcshapiro.test(D)
### T2 Hotelling Test 
# H0: delta == delta.0 vs H1: delta != delta.0
# with delta.0=c(0,0)
n <- dim(D)[1]  
p <- dim(D)[2]   
D.mean   <- colMeans(D)
D.cov    <- cov(D)
D.invcov <- solve(D.cov)
alpha   <- .05
delta.0 <- c(0,0)

D.T2 <- n * (D.mean-delta.0) %*% D.invcov %*% (D.mean-delta.0)
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
D.T2 < cfr.fisher #rifiuto se FALSE
P <- 1-pf(D.T2*(n-p)/(p*(n-1)), p, n-p)#p-value

### Simultanouse T2 intervals
IC.T2.DBOD <- c( D.mean[1]-sqrt(cfr.fisher*D.cov[1,1]/n) , D.mean[1], D.mean[1]+sqrt(cfr.fisher*D.cov[1,1]/n) )
IC.T2.DSS  <- c( D.mean[2]-sqrt(cfr.fisher*D.cov[2,2]/n) , D.mean[2], D.mean[2]+sqrt(cfr.fisher*D.cov[2,2]/n) )

### Bonferroni intervals
k <- p  # 2
cfr.t <- qt(1-alpha/(2*k),n-1)

IC.BF.DBOD <- c( D.mean[1]-cfr.t*sqrt(D.cov[1,1]/n) , D.mean[1], D.mean[1]+cfr.t*sqrt(D.cov[1,1]/n) )
IC.BF.DSS  <- c( D.mean[2]-cfr.t*sqrt(D.cov[2,2]/n) , D.mean[2], D.mean[2]+cfr.t*sqrt(D.cov[2,2]/n) )
#vettoriale:
IC.BF <- cbind( D.mean-cfr.t*sqrt(diag(D.cov)/n) , D.mean, D.mean+cfr.t*sqrt(diag(D.cov)/n) )
dimnames(IC.BF)[[2]] <- c ('inf','centr','sup')


# Test: H0.i: mu.i1 == mu.i2  vs H1.i: mu.i1 != mu.i2
z.i <- (x.mean1-x.mean2)/sqrt(x.var*(1/n1+1/n2))
p.i <- ifelse(z.i<0, 2*pnorm(z.i),2*(1-pnorm(z.i)))
which(p.i<.01)
# Bonferoni test: controlla prob di rifiutare H0 quando vera
k <- 520
which(p.i*k<.01)  
# Benjamini-Hockberg (control the false discovery rate)  
p.BH <- p.adjust(p.i, method='BH')



###LAB-6: Test for repeated measures ----
n <- dim(pressure)[1] #50
q <- dim(pressure)[2] #4

mcshapiro.test(pressure)
matplot(t(pressure), type='l')

M <- sapply(pressure,mean)
S <- cov(pressure)
# build one of the possible contrast matrices (q-1)xq
C <- matrix(c(-1, 1, 0, 0,
              -1, 0, 1, 0,
              -1, 0, 0, 1), 3, 4, byrow=T)
# Test: H0: C%*%mu == 0 vs H1: C%*%mu != 0
alpha   <- .05
delta.0 <- c(0,0,0)
Md <- C %*% M 
Sd <- C %*% S %*% t(C)
Sdinv <- solve(Sd)

T2 <- n * t( Md - delta.0 ) %*% Sdinv %*% ( Md - delta.0 )
cfr.fisher <- ((q-1)*(n-1)/(n-(q-1)))*qf(1-alpha,(q-1),n-(q-1)) 
T2 < cfr.fisher
P <- 1-pf(T2*(n-(q-1))/((q-1)*(n-1)),(q-1),n-(q-1)) #p-value

# Simultaneous T2 intervals
IC.T2 <- cbind( Md - sqrt(cfr.fisher*diag(Sd)/n) , Md, Md + sqrt(cfr.fisher*diag(Sd)/n) )
dimnames(IC.T2)[[2]] <- c('inf','sup')

# Bonferroni intervals 
k     <- q - 1   # number of increments (i.e., dim(C)[1])
cfr.t <- qt(1-alpha/(2*k),n-1)

IC.BF <- cbind( Md - cfr.t*sqrt(diag(Sd)/n) , Md, Md + cfr.t*sqrt(diag(Sd)/n) )
#Verifico per ogni intervallo se è contenuto lo 0:
# - se 0 è contenuto accetto H0
# - se 0 non contenuto, rifiuto H0
for (i in 1:q-1)
  print(paste('Reject H0 for a',i,': ', !(0>ICT2[i,1] & 0<ICT2[i,3]),sep=''))



###LAB-6: Test for two independent Gaussian populations ----
#data t1 & t2
n1 <- dim(t1)[1] # n1=3
n2 <- dim(t2)[1] # n2=4
p  <- dim(t1)[2] # p=2

t1.mean <- sapply(t1,mean)
t2.mean <- sapply(t2,mean)
t1.cov  <-  cov(t1)
t2.cov  <-  cov(t2)
Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2)
# Compare the matrices
list(S1=t1.cov, S2=t2.cov, Spooled=Sp)
# Test H0: mu1-mu2 == c(0,0)  vs  H1: mu1-mu2 != c(0,0)
alpha   <- .01
delta.0 <- c(0,0)
Spinv   <- solve(Sp)

T2 <- n1*n2/(n1+n2) * (t1.mean-t2.mean-delta.0) %*% Spinv %*% (t1.mean-t2.mean-delta.0)

cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
T2 < cfr.fisher # reject if FALSE at level alpha

P <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)

# Simultaneous T2 intervals
IC.T2.X1 <- c(t1.mean[1]-t2.mean[1]-sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)), t1.mean[1]-t2.mean[1]+sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)) )
IC.T2.X2 <- c(t1.mean[2]-t2.mean[2]-sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)), t1.mean[2]-t2.mean[2]+sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)) )
IC.T2 <- rbind(IC.T2.X1, IC.T2.X2)
dimnames(IC.T2)[[2]] <- c('inf','sup')                        
IC.T2

#Test for the mean & variance of two indipendent variance pop.
t.test(sardine[,1], sardine[,2], var.eq=T)#medie diverse
var.test(sardine[,1], sardine[,2]) #varianze uguale

#Bonferroni
k <- 2
cfr.t<-qt(1-alpha/(2*k),n1+n2-2)
IC.BN <- cbind(inf=t1.mean-t2.mean-cfr.t*sqrt(diag(Sp)*(1/n1+1/n2)), 
               sup=t1.mean-t2.mean+cfr.t*sqrt(diag(Sp)*(1/n1+1/n2)))
IC.BN



###LAB-7: One-way ANOVA  ----
#One way ANOVA (p=1,g=6)
x11()
plot(feed, weight, xlab='treat', ylab='weight', col='grey85', main='Dataset Chicken Weights')
### Model: weigth.ij = mu + tau.i + eps.ij; eps.ij~N(0,sigma^2)
### Test:
### H0: tau.1 = tau.2 = tau.3 = tau.4 = tau.5 = tau.6 = 0 vs H1: (H0)^c
### H0: The feed supplements don't have effect

#grafico media sotto h0 e sotto h1
x11()
par(mfrow=c(1,2))
barplot(rep(mean(weight),6), names.arg=levels(feed), ylim=c(0,max(weight)),
        las=2, col='grey85', main='Model under H0')
barplot(tapply(weight, feed, mean), names.arg=levels(feed), ylim=c(0,max(weight)),
        las=2, col=rainbow(6),main='Model under H1')

n       <- length(feed)      # total number of obs.
ng      <- table(feed)       # number of obs. in each group
treat   <- levels(feed)      # levels of the treatment
g       <- length(treat)     # number of levels (i.e., of groups)

### verify the assumptions:
# 1) normality (univariate) in each group (6 tests)
Ps <- c(shapiro.test(weight[ feed==treat[1] ])$p,
        shapiro.test(weight[ feed==treat[2] ])$p) #devono essere tutti e 6 i treat

# 2) same covariance structure (= same sigma^2)
Var <- c(var(weight[ feed==treat[1] ]),
         var(weight[ feed==treat[2] ]))

# test of homogeneity of variances
# H0: sigma.1 = sigma.2 = sigma.3 = sigma.4 = sigma.5 = sigma.6 VS H1: there exist i,j s.t. sigma.i!=sigma.j
bartlett.test(weight, feed) 

#ANOVA
fit <- aov(weight ~ feed)
summary(fit) #se rifiuto i trattamenti hanno effetto

#Bonferroni
k <- g*(g-1)/2 #n. of comparison
alpha= 0.05

Mediag  <- tapply(weight, feed, mean)#variabile da stimare,gruppo
SSres <- sum(residuals(fit)^2)
S <- SSres/(n-g)
ng <- c(sum(area=="Cancun"),sum(area=="Guanajato"),sum(area=="MexicoCity"))

ICrange=NULL
for(i in 1:(g-1)) {
  for(j in (i+1):g) {
    print(paste(levels(treat)[i],"-",levels(treat)[j]))        
    print(as.numeric(c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                       Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
    ICrange=rbind(ICrange,as.numeric(c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                                       Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
  }}

#grafico
x11(width = 14, height = 7)
par(mfrow=c(1,2))
plot(feed, weight, xlab='treat', ylab='weight', col = rainbow(6), las=2)

#plotto le differenze per ogni coppia
h <- 1
plot(c(1,g*(g-1)/2),range(ICrange), pch='',xlab='pairs treat', ylab='Conf. Int. tau weight')
for(i in 1:(g-1)) {
  for(j in (i+1):g) {
    ind <- (i-1)*g-i*(i-1)/2+(j-i)
    lines (c(h,h), c(ICrange[ind,1],ICrange[ind,2]), col='grey55'); 
    points(h, Mediag[i]-Mediag[j], pch=16, col='grey55'); 
    points(h, ICrange[ind,1], col=rainbow(6)[j], pch=16); 
    points(h, ICrange[ind,2], col=rainbow(6)[i], pch=16); 
    h <- h+1
  }}
abline(h=0)
#chi non contiene lo zero è significativo

##BONFERRONI
n <- dim(museo)[1]
g <- 3
k <- g*(g-1)/2+1
S <- sum(residuals(fit3)^2)/(n-g)

alpha<- .1

Mg  <- tapply(museo[,1], tipo, mean) 

label <- levels(factor(tipo))
n1 <- length(museo[tipo==label[1],1])
n2 <- length(museo[tipo==label[2],1])
n3 <- length(museo[tipo==label[3],1])
t <- qt(1-alpha/(2*k),n-g)

# Conf int for the means
ICB1<-data.frame(inf=Mg[1]-sqrt(S*(1/n1))*t,centr=Mg[1],sup=Mg[1]+sqrt(S/n1)*t)
ICB2<-data.frame(inf=Mg[2]-sqrt(S*(1/n2))*t,centr=Mg[2],sup=Mg[2]+sqrt(S/n2)*t)
ICB3<-data.frame(inf=Mg[3]-sqrt(S*(1/n3))*t,centr=Mg[3],sup=Mg[3]+sqrt(S/n3)*t)
ICB<-data.frame(rbind(ICB1,ICB2,ICB3))
ICB

# Conf int for variances
chi_u <- qchisq(alpha/(2*k),n-g)
chi_l <- qchisq(1-alpha/(2*k),n-g)
ICBV <- data.frame(inf=(n-g)*S/chi_l,centr=S,sup=(n-g)*S/chi_u)
ICBV


###LAB-7: One-way MANOVA ----
#(p=4, g=3), p=variables, g=groups
iris4 <-iris[,1:4]#tolgo la colonna con il treatment
### Model: X.ij = mu + tau.i + eps.ij; eps.ij~N_p(0,Sigma), X.ij, mu, tau.i in R^4
### Test: H0: tau.1 = tau.2 = tau.3  = (0,0,0)' VS H1: (H0)^c
i1 <- which(species.name=='setosa')
n1 <- length(i1)
n2 <- length(i2)
n3 <- length(i3)
n  <- n1+n2+n3
g  <- length(levels(species.name))
p  <- 4
group.name <-dataset[,1]
### Verify the assumptions:
# 1)  normality (multivariate) in each group (3 tests)
Ps <- NULL
for(i in 1:g)
  Ps <- c(Ps, mcshapiro.test(iris[get(paste('i',i, sep='')),1:4])$p) 
Ps #un p-value per ogni levello del gruppo (g)

# 2) same covariance structure (= same covariance matrix Sigma)
S  <-  cov(iris4)
S1 <-  cov(iris4[i1,]) #primo treat
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

fit <- manova(as.matrix(iris4) ~ species.name)
summary.manova(fit,test="Wilks") #if reject there is stat.evid to state that factor "species" has an effect on the mean feature on the flower
#Who's the responsible? vediamolo via ANOVA e Bonf.
# Exact tests for p<=2 or g<=3 already implemented in R

#Via ANOVA: for each of the p=4 variables we perform an ANOVA test to verify if the membership to a group has influence
# on the mean of the variable:
summary.aov(fit)

### Via Bonferroni
alpha <- 0.05
k <- p*g*(g-1)/2 #n. of comparison
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



###LAB-7: Two-way ANOVA -----
##### (p=1, g=2, b=2)
### => we verify the assumptions on the reduced model (one-way ANOVA)
# 1) normality (univariate) in each group (2 tests) 
Ps <- c(shapiro.test(durata[ AR=='aero_centro' & FF=='festivo' ])$p,
        shapiro.test(durata[ AR=='aero_centro' & FF=='feriale' ])$p,
        shapiro.test(durata[ AR=='centro_aero' & FF=='festivo' ])$p,
        shapiro.test(durata[ AR=='centro_aero' & FF=='feriale' ])$p)
Ps
#(durata--> variabili su cui faccio anoca)
# 2) homogeneity of variances
bartlett.test(list(durata[ AR=='aero_centro' & FF=='festivo' ],
                   durata[ AR=='aero_centro' & FF=='feriale' ],
                   durata[ AR=='centro_aero' & FF=='festivo' ],
                   durata[ AR=='centro_aero' & FF=='feriale' ]))

#Oppure --> metto as.factor tutto:
type.food <- as.factor(type.food)
area <- as.factor(area)
bartlett.test(price,type.food:area)
#price quella che voglio stimare

distr       <- factor(c('Esso','Esso','Esso','Esso','Shell','Shell','Shell','Shell'))
benz        <- factor(c('95','95','98','98','95','95','98','98'))
distr_benz  <- factor(c('Esso95','Esso95','Esso98','Esso98','Shell95','Shell95','Shell98','Shell98'))

g <- length(levels(distr)) #2
b <- length(levels(benz)) #2
n <- length(km)/(g*b) #2

#grafico con tutti i modelli
M           <- mean(km)
Mdistr      <- tapply(km,      distr, mean)
Mbenz       <- tapply(km,       benz, mean)
Mdistr_benz <- tapply(km, distr_benz, mean)

x11()
par(mfrow=c(2,3),las=2)
barplot(rep(M,4), names.arg=levels(distr_benz), ylim=c(0,24), main='No factor')
barplot(rep(Mdistr,each=2), names.arg=levels(distr_benz), ylim=c(0,24), 
        col=rep(c('blue','red'),each=2), main='Only Fact. Stat.')
barplot(rep(Mbenz,times=2), names.arg=levels(distr_benz), ylim=c(0,24),
        col=rep(c('darkgreen','orange'),times=2), main='Only Fact. Gas')
barplot(c(Mdistr[1]+Mbenz[1]-M, Mdistr[1]+Mbenz[2]-M, Mdistr[2]+Mbenz[1]-M, 
          Mdistr[2]+Mbenz[2]-M), names.arg=levels(distr_benz), ylim=c(0,24), 
        col=rep(c('darkgreen','orange'),times=2), density=rep(10,4), angle=135, 
        main='Additive model Stat.+Gas')
barplot(c(Mdistr[1]+Mbenz[1]-M, Mdistr[1]+Mbenz[2]-M, Mdistr[2]+Mbenz[1]-M, 
          Mdistr[2]+Mbenz[2]-M), names.arg=levels(distr_benz), ylim=c(0,24), 
        col=rep(c('blue','red'),each=2), density=rep(10,4), add=T)
barplot(Mdistr_benz, names.arg=levels(distr_benz), ylim=c(0,24), 
        col=rainbow(5)[2:5], main='Model with Interact. Stat.+Gas.')
plot(distr_benz, km, col=rainbow(5)[2:5], ylim=c(0,24),xlab='')

### Model with interaction (complete model): 
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N(0,sigma^2), i=1,2 (effect station), j=1,2 (effect gasoline)
#tau : effect stat, beta: effect gas, gamma: interaction between stat and gas
fit <- aov(km ~ distr + benz + distr:benz)
summary.aov(fit)
#rimuovo uno alla volta le cose non significative --> tolgo interazione gamma
fit <- aov(km ~ distr + benz)
summary.aov(fit)
# Remark: by removing the interaction, the residual degrees of freedom increase! 

### Test:
### 1) H0: gamma.11 = gamma.12 = gamma.21 = gamma.22 = 0    vs   H1: (H0)^c
###    i.e.,
###    H0: There is no significant interaction between the factors station
###        and gasoline in terms of performances
###    H1: There exists a significant interaction between the factors station 
###        and gasoline in terms of performances
###
### 2) H0: tau.1 = tau.2 = 0    vs   H1: (H0)^c
###    i.e.,
###    H0: The effect "gas station" doesn't significantly influence performances 
###    H1: The effect "gas station" significantly influences performances
###
### 3) H0: beta.1 = beta.2 = 0    vs   H1: (H0)^c
###    i.e.,
###    H0: The effect "gasoline" doesn't significantly influence performances
###    H1: The effect "gasoline" significantly influences performances

### Example: global test for the significance of the two treatments 
###          (model without interaction)
SSdistr <- sum(n*b*(Mdistr - M)^2)              # or from the summary: 1.53    
SSbenz  <- sum(n*g*(Mbenz  - M)^2)              # or from the summary: 66.70
SSres   <- sum((km - M)^2) - (SSdistr+SSbenz)   # or from the summary: 16.37

Ftot      <- ( (SSdistr + SSbenz) / ((g-1)+(b-1)))/(SSres / (n*g*b-g-b+1))
Ptot      <- 1 - pf(Ftot, (g-1)+(b-1), n*g*b-g-b+1) # attention to the dof!

#Riduco ad ANOVA --> per vedere gli effetti singoli
fit.aov1 <- aov(km ~ benz)
summary.aov(fit.aov1)

#Which is the best type of gasoline?
SSres <- sum(residuals(fit.aov1)^2)

### Interval at 90% for the differences (reduced additive model)
### [b=2, thus one interval only]
IC <- c(diff(Mbenz) - qt(0.95, (n*g-1)*b) * sqrt(SSres/((n*g-1)*b) *(1/(n*g) + 1/(n*g))), 
        diff(Mbenz) + qt(0.95, (n*g-1)*b) * sqrt(SSres/((n*g-1)*b) *(1/(n*g) + 1/(n*g))))
names(IC) <- c('Inf', 'Sup')
IC

# Estimate variances
W <- sum(fit$residuals^2)  # SS_res
var <- W/(g*b*n-g-b+1)     # SS_res/gdl(res)
var

# Estimate the great mean mu:
m <- mean(km)

# Estimate tau.i, beta.j:
tauAC  <- mean(euros[euros$AR=='aero_centro',1]) - m  # tau.1
tauCA  <- mean(euros[euros$AR=='centro_aero',1]) - m  # tau.2

betaFest <- mean(euros[euros$FF=='festivo',1]) - m  # beta.1
betaFer  <- mean(euros[euros$FF=='feriale',1]) - m  # beta.2

# Point-wise estimates of mean duration of travels
# (model without interaction!)
mAC_Fest <- m + tauAC + betaFest
mAC_Fer  <- m + tauAC + betaFer
mCA_Fest <- m + tauCA + betaFest
mCA_Fer  <- m + tauCA + betaFer


#BONFERRONI PER IL MODELLO RIDOTTO: means adn variance
N <- dim(doughnuts)[1]
g <- length(levels(tipo))
DF <- N-g

alpha <- .05
k <- g+1

qT <- qt(1-alpha/(2*k), DF)
qCinf <- qchisq(1 - alpha / (2*k), DF)
qCsup <- qchisq(alpha / (2*k), DF)

Spooled <- (t(fit.c3$res) %*% fit.c3$res)/DF   
Spooled

m1 <- mean(doughnuts[which(tipo==levels(factor(tipo))[1]),1])
m2 <- mean(doughnuts[which(tipo==levels(factor(tipo))[2]),1])
m3 <- mean(doughnuts[which(tipo==levels(factor(tipo))[3]),1])
medie <- c(m1,m2,m3)

ng <- c(length(which(tipo==levels(factor(tipo))[1])),
        length(which(tipo==levels(factor(tipo))[2])),
        length(which(tipo==levels(factor(tipo))[3])))

BF    <- rbind(cbind(inf=medie - sqrt(c(Spooled) / ng) * qT,
                     sup=medie + sqrt(c(Spooled) / ng) * qT),
               c(inf=Spooled * DF / qCinf,
                 sup=Spooled * DF / qCsup))
BF




###LAB-7: Two-way MANOVA ----
##### (p=3, g=2, b=2)
Ex   <- factor(plastic$Ex, labels=c('L','H')) # Treat.1
Ad   <- factor(plastic$Ad, labels=c('L','H')) # Treat.2

ExAd <- Ex
levels(ExAd) <- c('LL','LH','HL','HH')
ExAd[Ex=='L' & Ad=='L'] <- 'LL'
ExAd[Ex=='L' & Ad=='H'] <- 'LH'
ExAd[Ex=='H' & Ad=='L'] <- 'HL'
ExAd[Ex=='H' & Ad=='H'] <- 'HH'
#plastic3: tengo solo le colonne p
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

### Model with interaction (complete model): 
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N_p(0,Sigma), [p=3] i=1,2 (effect Extrusion), j=1,2 (effect Additive),
### Verify the assumptions (although we only have 5 data in each group!)
# 1) normality (multivariate) in each group (4 test)
Ps <- c(mcshapiro.test(plastic3[ ExAd==levels(ExAd)[1], ])$p,
        mcshapiro.test(plastic3[ ExAd==levels(ExAd)[2], ])$p,
        mcshapiro.test(plastic3[ ExAd==levels(ExAd)[3], ])$p,
        mcshapiro.test(plastic3[ ExAd==levels(ExAd)[4], ])$p)

# 2) homogeneity of the covariance (qualitatively)
#ATTENZIONE PRENDO IL DATSET CON I TRE PARAMETRI
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
#pattern different --> however we used assumption of homogeneous

### Model with interaction (complete model): 
fit <- manova( as.matrix(plastic3) ~ Ex + Ad + Ex:Ad)
summary.manova(fit, test="Wilks")

### Model without interaction (additive model): 
fit2<- manova( as.matrix(plastic3) ~ Ex + Ad)
summary.manova(fit2, test="Wilks")

# ANOVA on the components (we look at the 3 axes-directions in R^3 separately)
#Effect are present in the three components?
summary.aov(fit2)
#Response Tr : si Ex, si Ad --> effect of Ex and Ad in this direction
#Response Gl : si Ex, si Ad (10%)
#Response Op : no Ex, no Ad --> effect of treatment not significant

# Bonferroni: how much the extrusion changes the three resistence?
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
infEx <- mExH-mExL - qT * sqrt( diag(W)/(g*b*n-g-b+1) * (1/(n*g)+1/(n*g)) )
supEx <- mExH-mExL + qT * sqrt( diag(W)/(g*b*n-g-b+1) * (1/(n*g)+1/(n*g)) )

mAdL  <- sapply(plastic3[Ad=='L',],mean)
mAdH  <- sapply(plastic3[Ad=='H',],mean)
infAd <- mAdH-mAdL - qT * sqrt( diag(W)/(g*b*n-g-b+1) * (1/(n*g)+1/(n*g)) )
supAd <- mAdH-mAdL + qT * sqrt( diag(W)/(g*b*n-g-b+1) * (1/(n*g)+1/(n*g)) )

IC2   <- list(ExH_ExL=cbind(infEx, supEx), AdH_AdL=cbind(infAd, supAd))
IC2






###LAB-8: Linear and Quadratic Discriminant Analysis ----
library(MASS)
#Separare il dataset
si_snow <- snow[which(snow$Snow.G.1=="snow"),1:2]
no_snow <- snow[which(snow$Snow.G.1=="no-snow"),1:2]
group <-as.factor(snow[,3])
### LDA (univariate) 
# Assumptions:
# 1) if L=i, X.i ~ N(mu.i, sigma.i^2), i=A,B
# 2) sigma.A=sigma.B
# 3) c(A|B)=c(B|A) (equal misclassification costs)
# 1) normality (univariate) within the groups
shapiro.test(cyto[A,1]) #or mc.shapiro
shapiro.test(cyto[B,1])
# 2) equal variance (univariate)
var.test(cyto[A,1],cyto[B,1]) #or bartlett.test(data_senzaclassificatore, gruppi)

nA <- length(A)
nB <- length(B)
n  <- nA + nB

# Prior probabilities (estimated from the data, no prior knowledge)
PA <- nA / n
PB <- nB / n

cyto.lda <- lda(dataset,group)
cyto.lda

# misclassification costs
# misclassification costs
c.vf <- 10 #banconota falsa riconosciuta come vera
c.fv <- 0.05

pf <- 0.001
pt <- 1-0.001

# Prior modified to account for the misclassification costs
prior.c <- c(true=pt*c.fv/(c.vf*pf+c.fv*pt), false=pf*c.vf/(c.vf*pf+c.fv*pt))
prior.c


# posterior probability and classification for x=0
x <- data.frame(Infg = 0) #prova <- as.data.frame(cbind(quantita=200,temperatura=-4)) 
# The command predict() returns a list containing (see the help of predict.lda):
# - the class associated with the highest posterior probability 
predict(cyto.lda, x)$class
# - the posterior probabilities for the classes
predict(cyto.lda, x)$posterior
# - in lda: the coordinates of the canonical analysis of Fisher
#           (Fisher's discriminant scores)
predict(cyto.lda, x)$x

# set prior probabilities
cyto.lda.1 <- lda(group ~ Infg, prior=c(0.95,0.05))

#k-nearest neighbor classifier
library(class)
cyto.knn <- knn(train = Infg, test = x, cl = group, k = 3, prob=T)
cyto.knn.class <- (cyto.knn == 'B')+0 
cyto.knn.B <- ifelse(cyto.knn.class==1, 
                     attributes(cyto.knn)$prob, 
                     1 - attributes(cyto.knn)$prob)

x11()
plot(x[,1], cyto.LDA.B, type='l', col='red', lty=2, xlab='x', ylab='estimated posterior')
points(x[,1], cyto.knn.B, type='l', col='black', lty=1)
abline(h = 0.5)
legend(-10, 0.75, legend=c('LDA','knn'), lty=c(2,1), col=c('red','black'))

### k-nearest neighbor classifier in CV
library(class)
set.seed(123)
AER_CV<-numeric(21)
for(k in  10:30){
  a<-knn.cv(train = debris[,-3], cl = group, k = k, prob = FALSE, use.all = TRUE)
  misc<-table(class.true=group, class.assigned=a)
  G=2
  for(g in 1:G)
    AER_CV[k-9] <- AER_CV[k-9]+ sum(misc[g,-g])/sum(misc[g,]) 
  
}
min(AER_CV)#0.2013575 --> K=22
#Plot knn best region:
x  <- seq(min(debris[,1]), max(debris[,1]), length=300)
y  <- seq(min(debris[,2]), max(debris[,2]), length=300)
xy <- expand.grid(x=x, y=y)
data.knn <- knn(train = debris[,-3], test = xy, cl =group, k =22  )
z  <- as.numeric(data.knn)
x11()
plot(data_red,pch=20)
points(high,col='red',pch=20)
points(low,col='green',pch=20)
contour(x, y, matrix(z, 300), levels=c(1.5, 2.5), drawlabels=F, add=T)

#prediction
knn.best<-knn(train=data[,-3],test=data.frame(x=1,y=-4),k=20,cl=data[,3])
knn.best

###LDA cross-validation
LdaCV.iris <- lda(iris2, species.name, CV=TRUE)  # CV cross validation
LdaCV.iris$class
species.name
table(class.true=species.name, class.assignedCV=LdaCV.iris$class)

errorsCV <- (LdaCV.iris$class != species.name)
errorsCV
sum(errorsCV)

AERCV   <- sum(errorsCV)/length(species.name)
AERCV

###QDA --> var/cov different in the group
# Assumptions:
# 1) if L=i, X.i ~ N(mu.i, sigma.i^2), i=A,B
qda.iris <- qda(iris2, species.name)
qda.iris
Qda.iris <- predict(qda.iris, iris2)

#Grafico classification region
fv <- as.factor(rep(c('vero','falso'), c(100,100)))
biglietti <- rbind(true, false)
x11()
plot(banknotes[,1:2], main='Banknotes', xlab='V1', ylab='V2', pch=20)
points(false, col='red', pch=20)
points(true, col='blue', pch=20)
legend('bottomleft', legend=levels(fv), fill=c('red','blue'), cex=.7)

points(qda1$means, pch=4,col=c('red','blue') , lwd=2, cex=1.5)

x  <- seq(min(banknotes[,1]), max(banknotes[,1]), length=200) #lenght= NUMERO OSSS
y  <- seq(min(banknotes[,2]), max(banknotes[,2]), length=200)
xy <- expand.grid(V1=x, V2=y)#V1 e V2 sostituisco con i veri nomi delle colonne

#se ho 3 classi
# z  <- predict(lda.iris, xy)$post  # these are P_i*f_i(x,y)  
# z1 <- z[,1] - pmax(z[,2], z[,3])  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
# z2 <- z[,2] - pmax(z[,1], z[,3])  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    
# z3 <- z[,3] - pmax(z[,1], z[,2])  # P_3*f_3(x,y)-max{P_j*f_j(x,y)}


z  <- predict(qda1, xy)$post  
z1 <- z[,1] - z[,2] 
z2 <- z[,2] - z[,1]  

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

# compute the APER
Qda_pred <- predict(qda1, iris2)
tab <-table(class.true=group, class.assigned=Qda_pred$class)
errorsq <- (Qda_pred$class != group)
APERq   <- sum(errorsq)/length(group)
APERq

Qda_pred <- predict(qda1, iris2)
misc <- table(class.true=group, class.assigned=Qda_pred$class)
G <- 3
prior=c(,)
APER <- 0
for(g in 1:G)
  APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior[g]
#Expected economic loss
(tab[1,2]*c.vf+tab[2,1]*c.fv)/60

# Compute the estimate of the AER by leave-out-out cross-validation 
QdaCV.iris <- qda(iris2, species.name, CV=T)
QdaCV.iris$class
species.name
table(class.true=species.name, class.assignedCV=QdaCV.iris$class)

errorsqCV <- (QdaCV.iris$class != species.name)
errorsqCV

AERqCV   <- sum(errorsqCV)/length(species.name)
AERqCV

#Trivial: vedo qual'è il gruppo con prior maggiore, il trivial=(gruppo_minore)/totale
### FISHER DISCRIMINANT ANALYSIS
#Assumptions: homogeneity of the covariance structure
#vedi riga 580

###LAB-8: Support Vector Machines ----
library(e1071)
#vedi 823

dat <- data.frame(x=data[,-3], y=group)
# To set the parameter C we can use the function tune(), which is based on cross-validation (10-fold)
set.seed (1)
tune.out <- tune(svm,y~.,data=dat2 ,kernel = 'linear',
                 ranges =list(cost=c(0.001 , 0.01, 0.1, 1,5,10,100) ))
summary(tune.out)
#the best cost is 4 because give the smallest error

bestmod <- tune.out$best.model
summary(bestmod)

plot(bestmod , dat, col =c('salmon', 'light blue'), pch=19, asp=1)

# Fit a Support Vector Machine (kernel = "radial") given a cost C
svmfit <- svm(y~., data=dat [train ,], kernel ='radial', gamma =1, cost =1)
summary(svmfit)

# Plot the SVM
x11()
plot(svmfit , dat, col =c('salmon', 'light blue'), pch=19, asp=1)

# Prediction for a new observation (command predict())
xtest <- matrix(rnorm (20*2) , ncol =2)
ytest <- sample(c(-1,1) , 20, rep=TRUE)
xtest[ytest ==1 ,] <- xtest[ytest ==1,] + 1
testdat <- data.frame(x=xtest , y=as.factor (ytest))

plot(xtest, col =ifelse(ytest==1, 'light blue', 'salmon'), 
     pch=19, xlab='x1', ylab='x2', asp=1)

ypred <- predict(bestmod,testdat)
table(true.label=testdat$y, assigned.label =ypred )



###LAB-9: Hierarchical clustering ----
library(mvtnorm)
library(rgl)
library(car)
#Dissimilarity matrix
iris.e <- dist(dataset, method='euclidean')
image(1:150,1:150,as.matrix(iris.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j')
#Osservo se i gruppi sono già visibilmente separati
plot(dataset)
iris.m <- dist(iris4, method='manhattan')
iris.c <- dist(iris4, method='canberra')

#clustering
iris.es <- hclust(iris.e, method='single')
iris.ea <- hclust(iris.e, method='average')
iris.ec <- hclust(iris.e, method='complete')
iris.ew <- hclust(iris.e, method='ward.D2')
#Dendogram (k = :  cluster)
x11()
par(mfrow=c(1,3))
plot(iris.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.es, k=2)
plot(iris.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.ec, k=2)
plot(iris.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.ea, k=2)

#taglia il dendogram
# Fix k=2 clusters:
cluster.ec <- cutree(iris.ec, k=2) # euclidean-complete
cluster.es <- cutree(iris.es, k=2) # euclidean-single
cluster.ea <- cutree(iris.ea, k=2) # euclidean-average

group1 <- run[cluster.ec=="1",]
group2 <- run[cluster.ec=="2",]  

n1 <- dim(group1)[1]
n2 <- dim(group2)[1]
n3 <- dim(group3)[1]

# interpret the clusters (misc=sample(150))
table(label.true = species.name[misc], label.cluster = cluster.es)
table(label.true = species.name[misc], label.cluster = cluster.ec)
table(label.true = species.name[misc], label.cluster = cluster.ea)

#biavariate: grafico punti
x11()
par(mfrow=c(1,3))

plot(dataset, main = 'Single linkage', col=cluster, pch=16, asp=1)
legend('topright',legend=c('1','2','3'),col=c(1,2,3),pch=16)
# compute the cophenetic matrices 
coph.es <- cophenetic(iris.es)
coph.ec <- cophenetic(iris.ec)
coph.ea <- cophenetic(iris.ea)

# compare with dissimilarity matrix (Euclidean distance)
x11()
layout(rbind(c(0,1,0),c(2,3,4)))
image(as.matrix(iris.e), main='Euclidean', asp=1 )
image(as.matrix(coph.es), main='Single', asp=1 )
image(as.matrix(coph.ec), main='Complete', asp=1 )
image(as.matrix(coph.ea), main='Average', asp=1 )
#più la cofanatic distance è simile alla matrice della distanza iniziale --> più il clustering
#è coerente con la distanza originale --> algoritmo coerente
#Se c'è separazione tra i gruppi --> meglio SL, altrimenti è affetto da chain effect

# compute cophenetic coefficients
es <- cor(iris.e, coph.es)
ec <- cor(iris.e, coph.ec)
ea <- cor(iris.e, coph.ea)

c("Eucl-Single"=es,"Eucl-Compl."=ec,"Eucl-Ave."=ea)#più alto meglio



#confident region + ellipse
n <- dim(group)[1]
p <-2
x.mean <- sapply(group,mean)
x.cov <- cov(group)
alpha <- 0.01
cfr.chisq <- qchisq(1-alpha,p)
#test gaussianity
x11()
plot(group, asp = 1, col='gold', pch=19, xlim=c(-10,50))
ellipse(center=x.mean, shape=x.cov/n, radius=(sqrt(cfr.chisq)),col = 'black', lty = 2, center.pch = 4)




###LAB-9: K-means clustering ----
result.k <- kmeans(dataset, centers=2) # Centers: fixed number of clusters
names(result.k)
result.k$cluster      # labels of clusters
result.k$centers      # centers of the clusters
result.k$totss        # tot. sum of squares
result.k$withinss     # sum of squares within clusters
result.k$tot.withinss # sum(sum of squares nei cluster)
result.k$betweenss    # sum of squares between clusters
result.k$size         # dimention of the clusters

### How to choose k:
### 1) evaluate the variability between the groups with respect to 
#   the variability within the groups
b <- NULL
w <- NULL
for(k in 1:10){
  result.k <- kmeans(Q, k)
  w <- c(w, sum(result.k$wit))
  b <- c(b, result.k$bet)
}
x11()
matplot(1:10, w/(w+b), pch='', xlab='clusters', ylab='within/tot', main='Choice of k', ylim=c(0,1))
lines(1:10, w/(w+b), type='b', lwd=2)
#take the one with the stronger decreasing






###LAB-9: Multidimensional Scaling ----
#Given the distances (dissimilarities) among n statistical units, look for the
#k-dimensional representation (k small) of the n statistical units such that
#the distances (dissimilarities) among the representations of the n units
# are as close as possible to the original distances (dissimilarities) among the n units.

location <- cmdscale(eurodist, k=2)

# Compute the "stress": the higher it is, the worse the matching between original distances and their
# geometrical representation through MDS
Stressk <- NULL
for(k in 1:4)
{
  location.k <- cmdscale(eurodist, k)
  Stress <- (sum( (as.vector(eurodist) - as.vector(dist(location.k)))^2)  /
               sum( as.vector(location.k)^2))^(1/2)
  Stressk <- c(Stressk, Stress) 
}

x11()
plot(1:4,Stressk,xlab='k',ylab='Stress',lwd=2)




###LAB-10: Linear Model ----
### Multiple linear regression
library(car)
D <- ifelse(fly[,2]=='M', 1, 0)  # dummy
#tolgo la riga age dal modello: tpd <- data.frame(tper[,-3],D=D)
#Model:
#distance = beta_0 + beta_1 * speed + beta_2 * speed^2 + Eps (linear in the parameters!)
### Assumptions:
## 1) Parameter estimation: E(Eps) = 0  and  Var(Eps) = sigma^2 
fit <- lm(distance ~ speed1 + speed2)
summary(fit) 
sum(residuals(fit)^2)/fit$df  # estimate of sigma^2

plot(cars, xlab='Speed', ylab='Stopping distance', las=1, xlim=c(0,30), ylim=c(-5,130))
x <- seq(0,30,by=0.1)
b <- coef(fit)
lines(x, b[1]+b[2]*x+b[3]*x^2)

#Grafico punti dati - fittati
x11()
plot(new_T$T,new_T$L,col='blue')
points(new_T$T, fit_2$fitted.values, col='green')

## 2) Inference:  Eps ~ N(0, sigma^2)
par(mfrow=c(2,2))
plot(fit)
shapiro.test(residuals(fit))
#Test (Fisher):
# H0: (beta1, beta2) == (0, 0) vs H1: (beta1, beta2) != (0, 0)
linearHypothesis(fit, rbind(c(0,1,0), c(0,0,1)), c(0,0)) 
summary(fm)

# d) Perform a test of level 5% to verify if during weekends the mean amount 
# of sold focaccia increases of 60 kg. In case, update the parameters’ estimates.
linearHypothesis(fit2, rbind(c(0,0,1)), c(60))
#p-value=0.85 --> i can't reject H0 
#true statement
peso <- kg - 60*d #constrained
fit3 <- lm(peso ~ t,data=focaccia2)
summary(fit3)

p <- 2  # number of tested coefficients
r <- 2  # number of regressors

# Confidence region:
# center: point estimate
c(coefficients(fm)[2], coefficients(fm)[3])
# Direction of the axes?
eigen(vcov(fm)[2:3,2:3])$vectors
x11()
plot(coefficients(fm)[2], coefficients(fm)[3], xlim = c(-6,6), ylim = c(-6,6), asp=1, xlab='beta1', ylab='beta2')
ellipse(coefficients(fm)[2:3], vcov(fm)[2:3,2:3], sqrt(p*qf(1-0.05,p,n-(r+1))))
abline(v=0)
abline(h=0)
# Note: collinearity! -->elongate shape of confidence region

# Bonferroni intervals (level 95%)
Bf <- rbind(
  beta1=c(coefficients(fm)[2]-sqrt(vcov(fm)[2,2])*qt(1-0.05/(2*p), n-(r+1)),
          coefficients(fm)[2]+sqrt(vcov(fm)[2,2])*qt(1-0.05/(2*p), n-(r+1))),
  beta2=c(coefficients(fm)[3]-sqrt(vcov(fm)[3,3])*qt(1-0.05/(2*p), n-(r+1)),
          coefficients(fm)[3]+sqrt(vcov(fm)[3,3])*qt(1-0.05/(2*p), n-(r+1))))
# or (only for intervals on beta)
confint(fm, level= 1-0.05/p)[2:3,]  # Bonferroni correction!


C <- rbind(c(1,0,0), c(1,1,0), c(1,1,1))
n <- dim(goat)[1]

Bf <- rbind(
  c((C %*% coefficients(fit.red))[1] - sqrt((C %*% vcov(fit.red) %*% t(C))[1,1]) * qt(1 - 0.10/6, n-3),
    (C %*% coefficients(fit.red))[1] + sqrt((C %*% vcov(fit.red) %*% t(C))[1,1]) * qt(1 - 0.10/6, n-3)),
  c((C %*% coefficients(fit.red))[2] - sqrt((C %*% vcov(fit.red) %*% t(C))[2,2]) * qt(1 - 0.10/6, n-3),
    (C %*% coefficients(fit.red))[2] + sqrt((C %*% vcov(fit.red) %*% t(C))[2,2]) * qt(1 - 0.10/6, n-3)),
  c((C %*% coefficients(fit.red))[3] - sqrt((C %*% vcov(fit.red) %*% t(C))[3,3]) * qt(1 - 0.10/6, n-3),
    (C %*% coefficients(fit.red))[3] + sqrt((C %*% vcov(fit.red) %*% t(C))[3,3]) * qt(1 - 0.10/6, n-3))
)
Bf


#Variance of error epsilon
e <- residuals(fm)
S2 <- t(e)%*%e / (df.residual(fm))
S2


##### Confidence intervals for the mean & prediction (new obs)
##### Assumption: Eps ~ N(0, sigma^2)
Z0.new <- data.frame(speed1=10, speed2=10^2)
alpha <- .1
k <- 2
n <- dim(index)[1]
r <- 3 #regressori
# Conf. int. for the mean
Conf <- predict(fit, Z0.new, interval='confidence', level=1-alpha/k)  
Conf
# Pred. int. for a new obs
Pred <- predict(fit, Z0.new, interval='prediction', level=1-alpha/k)  
Pred


#Conf. int for the variance
e <- residuals(fitB)
ICBvar <- data.frame(L=t(e)%*%e/qchisq(1-alpha/(2*k),n-(r+1)),
                     U=t(e)%*%e/qchisq(alpha/(2*k),n-(r+1)))
ICBvar

## Verify assumptions: homoschedasticity & gaussanity
par(mfrow=c(2,2))
plot(fm)
#1.homoschedasticity: (Residual vs Fitted): satisfy but there are some outliers in the central part
#if there is a trend --> heteroschedasticity
#2. Gaussanity: (QQplot & test)
shapiro.test(residuals(fm))

#DIAGNOSTIC
plot(X, Y, main='Scatterplot brain weight vs body weight', lwd=2,
     xlab='Body weight', ylab='Brain weight',col='white',xlim=c(-1000,8000))
text(X, Y,dimnames(data)[[1]],cex=1)
abline(coef[1],coef[2], lwd=2,col='red')

##Outlier
#--> eliminate
#log trasformation

### Model 2 (dummy variable)
Q <- cbind(quakes[,1:2], depth=-quakes[,3]/100)
d <- dist(Q)
clusterw <- cutree(hclust(d, method='ward.D2'), 2)

Qd <- cbind(Q, dummy)
head(Qd)
dummy <- clusterw - 1 

# Model: ANCOVA MODEL
# depth = beta0       + beta1*lat       + beta2*long        +
#       + beta3*dummy + beta4*dummy*lat + beta5*dummy*long  + eps
# i.e.,
# depth = B0[g] + B1[g]*lat + B2[g]*long + eps
# with B0[g]=beta0       if the unit is in group s.t. dummy=0 (red)
#      B0[g]=beta0+beta3 if the unit is in group s.t. dummy=1 (green)
#      B1[g]=beta1       if the unit is in group s.t. dummy=0 (red)
#      B1[g]=beta1+beta4 if the unit is in group s.t. dummy=1 (green)
#      B2[g]=beta2       if the unit is in group s.t. dummy=0 (red)
#      B2[g]=beta2+beta5 if the unit is in group s.t. dummy=1 (green)

fitd <- lm(depth ~ lat + long + dummy + lat:dummy + long:dummy, data=Qd)
summary(fitd)

# test: are the two planes needed?
A <- rbind(c(0,0,0,1,0,0), c(0,0,0,0,1,0), c(0,0,0,0,0,1))
b <- c(0,0,0)
linearHypothesis(fitd, A, b)

#Predict interval
new_data <- data.frame(Year2 = c(11,11,11)^2, dFr=c(1,0,0), dGer=c(0,1,0))
IP <- predict(fit2, newdata=new_data, interval='prediction', level=1-0.05/3)
rownames(IP) <- c('Fr','Ger','It')
IP


#Creare dummy variable
data <-cbind(temperature= c(Temperature[,1],Temperature[,2],Temperature[,3]),
             Sin = sin(2*pi/12*c(1:12,1:12,1:12)),
             Cos = cos(2*pi/12*c(1:12,1:12,1:12)),
             Res = c(rep(0,12),rep(1,12),rep(0,12)),
             Mon = c(rep(0,12),rep(0,12),rep(1,12)))

#Stampo coeff
Beta <- rbind(
  Edmonton = coef(fit)[c(1,4,5)],
  Resolute = coef(fit)[c(1,4,5)] +  coef(fit)[c(2,6,7)],
  Montreal = coef(fit)[c(1,4,5)] +  coef(fit)[c(3,8,9)])
Beta

###LAB-11: Problem with collinearity ----
library(MASS)
library(car)
library(rgl)
library(glmnet)
#Large p-value in all regressor --> collinearity
vif(fm) #>5/10 collinearity

##Solution: PCA
speed.pc <- princomp(cbind(speed1,speed2), scores=TRUE)
summary(speed.pc)
speed.pc$load

# Explained variance
layout(matrix(c(2,3,1,3),2,byrow=T))
barplot(result.pc$sdev^2, las=2, main='Principal Components', ylab='Variances')
barplot(c(sd(X1),sd(X2),sd(X3),sd(X4))^2, las=2, main='Original variables', ylab='Variances')
plot(cumsum(result.pc$sdev^2)/sum(result.pc$sde^2), type='b', axes=F, xlab='number of components', ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.9, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:4,labels=1:4,las=2)


# Loadings
par(mar = c(1,4,0,2), mfrow = c(4,1))
for(i in 1:4)barplot(result.pc$load[,i], ylim = c(-1, 1))
#contrast between x2 and x4 in comp.1
#contrast between x1 and x3 in comp.2

sp1.pc <- speed.pc$scores[,1]
sp2.pc <- speed.pc$scores[,2]
# Model: y = b0 + b1*PC1+ b2*PC2 + eps, eps~N(0,sigma^2)
fm.pc <- lm(distance ~ sp1.pc + sp2.pc)

# We can re-write the model as:
# Model: 
# y= b0 + b1*      PC1                 + b2*      PC2                 + eps =
#  = b0 + b1*(e11*(X1-m1)+e21*(X2-m2)) + b2*(e12*(X1-m1)+e22*(X2-m2)) + eps =
#  = b0 - b1*e11*m1 - b2*e12*m1 - b1*e21*m2 - b2*e22*m2 + 
#                           + (b1*e11+b2*e12)*X1 + (b1*e21+b2*e22)*X2 + eps
# where e.ij are the loadings, i=1,2, j=1,2.
# => We can compute the coefficients of the model which used the original 
#    regressors
m1 <- mean(X1)
m2 <- mean(X2)
m3 <- mean(X3)
m4 <- mean(X4)
beta0 <- coefficients(fm.pc)[1] -
  coefficients(fm.pc)[2]*result.pc$load[1,1]*m1 -
  coefficients(fm.pc)[3]*result.pc$load[1,2]*m1 -
  coefficients(fm.pc)[2]*result.pc$load[2,1]*m2 -
  coefficients(fm.pc)[3]*result.pc$load[2,2]*m2 -
  coefficients(fm.pc)[2]*result.pc$load[3,1]*m3 -
  coefficients(fm.pc)[3]*result.pc$load[3,2]*m3 - 
  coefficients(fm.pc)[2]*result.pc$load[4,1]*m4 -
  coefficients(fm.pc)[3]*result.pc$load[4,2]*m4
beta1 <- coefficients(fm.pc)[2]*result.pc$load[1,1] +
  coefficients(fm.pc)[3]*result.pc$load[1,2]
beta2 <- coefficients(fm.pc)[2]*result.pc$load[2,1] +
  coefficients(fm.pc)[3]*result.pc$load[2,2]
beta3 <- coefficients(fm.pc)[2]*result.pc$load[3,1] +
  coefficients(fm.pc)[3]*result.pc$load[3,2]
beta4 <- coefficients(fm.pc)[2]*result.pc$load[4,1] +
  coefficients(fm.pc)[3]*result.pc$load[4,2]

c(beta0=as.numeric(beta0),beta1=as.numeric(beta1),beta2=as.numeric(beta2),beta3=as.numeric(beta3),beta4=as.numeric(beta4))
result$coefficients




##Solution: ridge
set.seed(27289292)#IMPORTANTE
# Choice of the optimal lambda, e.g., via cross-validation
lambda.c <- seq(0,10,0.01)
fit.ridge <- lm.ridge(distance ~ speed1 + speed2, lambda = lambda)
select(fit.ridge)

# or
lambda.opt <- lambda.c[which.min(fit.ridge$GCV)]
lambda.opt

##Solution: Lasso

x <- model.matrix(distance~speed1+speed2)[,-1] # Build the matrix of predictors

y <- distance # Build the vector of response

# Let's set a grid of candidate lambda's for the estimate
lambda.grid <- 10^seq(5,-3,length=100)
fit.lasso <- glmnet(x,y, lambda = lambda.grid) # default: alpha=1 -> lasso 
# [note: if alpha=0 -> ridge regression]

plot(fit.lasso,xvar='lambda',label=TRUE, col = rainbow(dim(x)[2]))
legend('topright', dimnames(x)[[2]], col =  rainbow(dim(x)[2]), lty=1, cex=1)

# Let's set lambda via cross validation
cv.lasso <- cv.glmnet(x,y,lambda=lambda.grid) # default: 10-fold CV

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso #0.006428073

plot(cv.lasso)
abline(v=log(bestlam.lasso), lty=1)

# Get the coefficients for the optimal lambda
coef.lasso <- predict(fit.lasso, s=bestlam.lasso, type = 'coefficients')
coef.lasso 

plot(rep(0, dim(x)[2]), coef(lm(y~x))[-1], col=rainbow(dim(x)[2]), pch=20, xlim=c(-1,3), ylim=c(-1,2), xlab='', ylab=expression(beta),
     axes=F)
points(rep(1, dim(x)[2]), coef.ridge[-1], col=rainbow(dim(x)[2]), pch=20)
points(rep(2, dim(x)[2]), coef.lasso[-1], col=rainbow(dim(x)[2]), pch=20)
abline(h=0, col='grey41', lty=1)
box()
axis(2)
axis(1, at=c(0,1,2), labels = c('LS', 'Ridge', 'Lasso'))


###  Subset Selection Methods
library(leaps)
regfit.full <- regsubsets(Salary~., data=Hitters, nvmax=19)
summary(regfit.full)
reg.summary <- summary(regfit.full)
names(reg.summary)

x11(height=7,width=14)
par(mfrow=c(1,3))
plot(reg.summary$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
plot(reg.summary$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
plot(reg.summary$rss,xlab="Number of Variables",ylab="RSS",type="b")

which.max(reg.summary$adjr2)
coef(regfit.full,11)

x11()
plot(regfit.full,scale="r2",main="Exhaustive search")

x11()
plot(regfit.full,scale="adjr2",main="Exhaustive search")

#K-fold cross-validation

k <- 10

set.seed(1)
folds <- sample(1:k,nrow(Hitters),replace=TRUE)
folds
table(folds)

# function that performs the prediction for regsubsets
predict.regsubsets <- function(object,newdata,id){
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form,newdata)
  coefi <- coef(object,id=id)
  xvars <- names(coefi)
  mat[,xvars]%*%coefi
}

cv.errors <- matrix(NA,k,19, dimnames=list(NULL, paste(1:19)))
for(j in 1:k){
  best.fit <- regsubsets(Salary~.,data=Hitters[folds!=j,],nvmax=19)
  for(i in 1:19){
    pred <- predict(best.fit,Hitters[folds==j,],id=i)
    cv.errors[j,i] <- mean( (Hitters$Salary[folds==j]-pred)^2 )
  }
}
cv.errors
root.mean.cv.errors <- sqrt(apply(cv.errors,2,mean)) # average over the columns
root.mean.cv.errors
plot(root.mean.cv.errors,type='b')

which.min(root.mean.cv.errors)
points(which.min(root.mean.cv.errors),root.mean.cv.errors[which.min(root.mean.cv.errors)], col='red',pch=19)

### Ridge and Lasso regression with glmnet
x <- model.matrix(Salary~.,Hitters)[,-1] # predictor matrix
y <- Hitters$Salary # response
grid <- 10^seq(10,-2,length=100) # grid of lambda

# Ridge regression
ridge.mod <- glmnet(x,y,alpha=0,lambda=grid)#Per LASSO metto alpha=1
plot(ridge.mod,xvar='lambda',label=TRUE)

set.seed(123)
cv.out <- cv.glmnet(x,y,alpha=0,nfold=3,lambda=grid) 

bestlam.ridge <- cv.out$lambda.min
bestlam.ridge
log(bestlam.ridge)

#Coeff
coef.ridge <- predict(ridge.mod, s=bestlam.ridge, type = 'coefficients')[1:20,]
coef.lasso <- predict(lasso.mod, s=bestlam.lasso, type = 'coefficients')[1:20,]
coef.ridge
coef.lasso

plot(rep(0, dim(x)[2]), coef(lm(y~x))[-1], col=rainbow(dim(x)[2]), pch=20, xlim=c(-1,3), ylim=c(-1,2), xlab='', ylab=expression(beta),
     axes=F)
points(rep(1, dim(x)[2]), coef.ridge[-1], col=rainbow(dim(x)[2]), pch=20)
points(rep(2, dim(x)[2]), coef.lasso[-1], col=rainbow(dim(x)[2]), pch=20)
abline(h=0, col='grey41', lty=1)
box()
axis(2)
axis(1, at=c(0,1,2), labels = c('LS', 'Ridge', 'Lasso'))


###LAB-11: Classification and regression trees ----
tree.boston <- tree(medv~.,Boston)
summary(tree.boston)

plot(tree.boston)
text(tree.boston,pretty=0)

cv.boston <- cv.tree(tree.boston)

plot(cv.boston$size,cv.boston$dev,type='b',xlab='size',ylab='deviance')

prune.boston <- prune.tree(tree.boston,best=4)

plot(prune.boston)
text(prune.boston,pretty=0)








###GEOSTAT ----
#15_06_2020
library(sp)           
library(lattice)      
library(geoR)         
library(gstat)        
## Functions for graphics 
v.f <- function(x, ...){100-cov.spatial(x, ...)}
v.f.est<-function(x,C0, ...){C0-cov.spatial(x, ...)}

coordinates(data) <- c('Long','Lat')
names(data) <- c('Bq','Long','Lat','D')
##    Exploration analysis
bubble(dataset,'zinc',do.log=TRUE,key.space='bottom')
hist(meuse$zinc, breaks=16, col="grey", main='Histogram of Zn', prob = TRUE, xlab = 'Zn') #asymmetric
hist(log(meuse$zinc), breaks=16, col="grey", main='Histogram of log(Zn)', prob = TRUE, xlab = 'log(Zn)')
xyplot(log(zinc) ~ sqrt(dist), as.data.frame(meuse))

##    Variogram Analysis
v <- variogram(log(zinc) ~ 1, meuse) #model stationary ~1 (just used the intercept)
plot(v, main = 'Sample Variogram',pch=19)#conotrollo se si stabilizza
#if the estimation go faster than the quadratic --> problem in the estimation (problem in assumption?)
plot(variogram(log(zinc) ~ 1, meuse, alpha = c(0, 45, 90, 135)),pch=19)

##    Variogram modeling
# vgm(sill, model, range, nugget)
# weighted least squares fitting a variogram model to the sample variogram,steps:
# 1) choose a suitable model (both spherical and exponential model have a linear behavior near the
#         origin but exponential model has a faster growth than the spherical one)
# 2) choose suitable initial values for partial sill, range & nugget
## 3) fit the model using one of the possible fitting criteria
v.fit <- fit.variogram(v, vgm(1, "Sph", 800, 1))#sill: salto, .., range: quando conv
plot(v, v.fit, pch = 19)
#se non mi da il range lo decido guardando il plot della semivariance
## maximum likelihood fitting of variogram models
fit.variogram.reml(log(zinc)~1, meuse, model=vgm(0.6, "Sph", 800, 0.06))

# Linear behavior near the origin, growth not very fast 
# Recall: both spherical and exponential model have a linear behavior near the
#         origin but exponential model has a faster growth than the spherical one


##    SPATIAL PREDICTION & KRIGING 
## Stationary Univariate Spatial Prediction (Ordinary Kriging)
## Prediction in a single new location 
s0.new=data.frame(x=179180, y=330100) # UTM coordinates
#s0.new=as.data.frame(matrix(c(78.59,35.34,0),1,3))
coordinates(s0.new)=c('x','y')
g.tr <- gstat(formula = log(zinc) ~ 1, data = meuse, model = v.fit) #gstat(g.obj, id, formula, data, model, set,...)

##  Make the ordinary kriging prediction with the function: 
# predict(obj, grid, BLUE=FALSE)
predict(g.tr, s0.new)
#var1.pred: prediction at my new location
#var1.var: uncertaintly (solo in ordinary crikking!!!)
predict(g.tr, s0.new, BLUE = TRUE) #estimate of the mean (trend component) under gls
# prediction over the entire grid
lz.ok <- predict(g.tr, meuse.grid, BLUE = FALSE)
spplot(lz.ok[,2])
spplot(lz.ok)

## Non-stationary Univariate Spatial Prediction (Universal Kriging)
# gstat(g.obj, id, formula, data, model, set,...)
meuse.gstat <- gstat(id = 'zinc', formula = log(zinc) ~ sqrt(dist),
                     data = meuse, nmax = 50, model=v.fit, set = list(gls=1))
v.gls<-variogram(meuse.gstat)
plot(v.gls)

v.gls.fit <- fit.variogram(v.gls, vgm(1, "Sph", 800, 1))
plot(v.gls, v.gls.fit, pch = 19)
# Update gstat object with variogram model
meuse.gstat <- gstat(id = 'zinc', formula = log(zinc) ~ sqrt(dist),
                     data = meuse, nmax = 50, model=v.gls.fit, set = list(gls=1))
## universal kriging:
s0.vec <- as.vector(slot(s0.new,'coords'))
s0.dist <- min(rowSums(scale(meuse.riv,s0.vec)^2)) 
s0.new <- as.data.frame(c(s0.new,s0.dist))
names(s0.new) <- c('x','y','dist')
coordinates(s0.new) <- c('x','y')
s0.new <- as(s0.new, 'SpatialPointsDataFrame')

predict(meuse.gstat, s0.new)
predict(meuse.gstat, s0.new, BLUE = TRUE) #estimate of x(s_0)'*beta (trend component) under gls
lz.uk <- predict(meuse.gstat, meuse.grid, BLUE=FALSE)# prediction over the entire grid
lz.uk.BLUE <- predict(meuse.gstat, meuse.grid, BLUE=TRUE)# estimate of the mean over the entire grid
spplot(lz.ok[,1], main = 'Ordinary Kriging, gstat')
spplot(lz.uk[,1], main = 'Universal Kriging, gstat')
spplot(lz.uk.BLUE[,1], main = 'Universal Kriging - drift , gstat')


# c) Estimate, via Generalized Least Squares, the parameter(s) a of the model chosen at point (a).
s0new <- data.frame(x=402476,y=4605558,distance=0)#metto distanza zero cosi trovo a0, le coord non contano
s1new <- data.frame(x=402476,y=4605558,distance=1)#per trovare a1 predico dove dist=1 e poi sottraggo a0
coordinates(s0new) <- c('x','y')
coordinates(s1new) <- c('x','y')

g.tr <- gstat(formula = speed ~ distance, data = mont, model = v.fit)
predict(g.tr, s0new, BLUE = TRUE)

a0 <- predict(g.tr, s0new, BLUE = TRUE)$var1.pred
a1 <- predict(g.tr, s1new, BLUE = TRUE)$var1.pred - a0 
a0
a1



###PERMUTATION: univariate (exam 12_09_2019)(Exam 21_11_2019) ----
# We sample n1 data from a first population X1 and n2 data froma second population  X2
# Test:
# H0:  X1 =^d  X2
# H1:  X1 !=^d X2
# Parameters:
n1 <- n2 <- 10
n <- n1 + n2

permutation <- sample(1:n) 
x_perm <- x_pooled[permutation] 
x1_perm <- x_perm[1:n1]
x2_perm <- x_perm[(n1+1):n]

# Test statistic: absolute difference between the two means
T0 <- abs(mean(x1) - mean(x2))

# Cardinality of the permutational space:
factorial(n)

# Number of distinct values of T*:
factorial(n)/(2*factorial(n1)*factorial(n2))

# Minimun achieveable p-value:
1/(factorial(n)/(2*factorial(n1)*factorial(n2)))
x_pooled <- c(x1,x2)
# loop to estimate the p-value
for(perm in 1:B){
  # permutation:
  permutation <- sample(1:n)
  x_perm <- x_pooled[permutation]
  x1_perm <- x_perm[1:n1]
  x2_perm <- x_perm[(n1+1):n]
  # test statistic:
  T_stat[perm] <- abs(mean(x1_perm) - mean(x2_perm))
}

# Permutational distribution of T
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat))
abline(v=T0,col=3,lwd=2)

# p-value
p_val <- sum(T_stat>=T0)/B
p_val

#Da esame 12 sett
set.seed(123)
p_val<-numeric(8)
#ipotesi nulla mu0=mu1 per gruppo j
#H0:mu1==mu2 vs H1:mu1>mu2
mod <- stress[1:10,]
no_mod <- stress[11:15,]
n1 <- 10
n2 <-5
n <- n1+n2
for(i in 1:8){
  T0 <- abs(median(mod[,i]) - median(no_mod[,i]))
  stat[i] <-T0
  T_stat <- numeric(10000)
  x_poled <- c(mod[,i],no_mod[,i])
  for(perm in 1:10000){
    permutation <- sample(1:n)
    x_perm <- x_poled[permutation]
    x1_perm <- x_perm[1:n1]
    x2_perm <- x_perm[(n1+1):n]
    # test statistic:
    T_stat[perm] <- abs(median(x1_perm) - median(x2_perm))
  }
  p_val[i] <-sum(T_stat>=T0)/10000
}
p_val # 0.1971 0.0730 0.0227 1.0000 0.0758 0.0661 0.1167 0.3378

T0 <-numeric(8)
for(i in 1:8){
  T0[i] <- abs(median(mod[,i]) - median(no_mod[,i]))
}  
T0

#False discovery rate
k<-8
alpha<-0.25
p<-sort(p_val)
vec<-(1:8)*alpha/k#sequenza di valori m=max{j:p(j)<=j/k*alpha} k valori di p_val
#in questo caso vorrà dire che c'è differenza , il secondo fornitore ha caratteristiche migliori 

# Bonferoni test: controlla prob di rifiutare H0 quando vera
p.Bf <- p.adjust(p_val, method='bonferroni')
which(p.Bf<.01)
###PERMUTATION: multivariate ----
##Two (independent) multivariate population test
t1 <- week[1:5,]
t2 <- week[6:7,]

#H0: non c'è differenza tra weekend e giorni lavorativi
# Computing a proper test statistic 
# (i.e., squared distance between the two sample mean vectors)
t1.mean <- colMeans(t1)
t2.mean <- colMeans(t2)

n1 <- dim(t1)[1]
n2 <- dim(t2)[1]
n  <- n1 + n2

T20 <- as.numeric((t1.mean-t2.mean) %*% (t1.mean-t2.mean))
T20

# Selecting a proper permutation strategy
# (i.e., data point permutations)

# number of possible data point permutations 
factorial(7)
# number of different values of the test statistic
choose(7,5)

# Estimating the permutational distribution under H0
B <- 100000
T2 <- numeric(B)

for(perm in 1:B){
  # Random permutation of indexes
  # When we apply permutations in a multivariate case, we keep the units together
  # i.e., we only permute the rows of the data matrix
  t_pooled <- rbind(t1,t2)
  permutation <- sample(n)
  t_perm <- t_pooled[permutation,]
  t1_perm <- t_perm[1:n1,]
  t2_perm <- t_perm[(n1+1):n,]
  
  # Evaluation of the test statistic on permuted data
  t1.mean_perm <- colMeans(t1_perm)
  t2.mean_perm <- colMeans(t2_perm)
  T2[perm]  <- (t1.mean_perm-t2.mean_perm) %*% (t1.mean_perm-t2.mean_perm) 
}

# plotting the permutational distribution under H0
hist(T2,xlim=range(c(T2,T20)),breaks=1000)

abline(v=T20,col=3,lwd=4)

plot(ecdf(T2))
abline(v=T20,col=3,lwd=4)

# p-value
p_val <- sum(T2>=T20)/B
p_val

##Center of simmetry of one multivariate population
##two paired multivariate population test
###PERMUTATION: ANOVA ---- 
###PERMUTATION: regression ----
###FDA: smoothing ----
library(fda)
#Control the derivative of the function, compute the central finite differences
rappincX1 <- (Xobs0[3:NT]-Xobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])
rappincX2 <- ((Xobs0[3:NT]-Xobs0[2:(NT-1)])/(abscissa[3:NT]-abscissa[2:(NT-1)])-(Xobs0[2:(NT-1)]-Xobs0[1:(NT-2)])/(abscissa[2:(NT-1)]-abscissa[1:(NT-2)]))*2/(abscissa[3:(NT)]-abscissa[1:(NT-2)])
par(mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(truecurve$Abscissa,truecurve$X0vera,type='l',col="orange",lwd=3)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(truecurve$Abscissa,truecurve$X1vera,type='l',col="orange",lwd=3)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(truecurve$Abscissa,truecurve$X2vera,type='l',col="orange",lwd=3)

### REGRESSION SPLINES 
m <- 5           # spline order 
degree <- m-1    # spline degree 
nbasis <- 9      #cross-validation for the best
#create the basis:
basis <- create.bspline.basis(rangeval=c(0,1), nbasis=nbasis, norder=m) 
plot(basis)

# Evaluate the basis on the grid of abscissa
basismat <- eval.basis(abscissa, basis)#,Lfdobj=k (k° derivata)
dim(basismat)
head(basismat)

# Fit via LS
lsfit(basismat, Xobs0, intercept=FALSE)$coef 
Xsp0 <- basismat %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef

par(mfrow=c(1,1))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
abline(v=basis$params)

# Approximate pointwise confidence intervals
# As in linear models, we can estimate the variance of x(t) as
# sigma^2*diag[phi*(phi'phi)^{-1}(phi)']
S <- basismat%*%solve(t(basismat)%*%basismat)%*%t(basismat) #projection operator 
sigmahat <- sqrt(sum((Xsp0-Xobs0)^2)/(NT-nbasis)) #estimate of sigma
lb <- Xsp0-qnorm(0.975)*sigmahat*sqrt(diag(S))
ub <- Xsp0+qnorm(0.975)*sigmahat*sqrt(diag(S))

plot(abscissa,Xsp0,type="l",col="blue",lwd=2,ylab="")
points(abscissa,lb,type="l",col="blue",lty="dashed")
points(abscissa,ub,type="l",col="blue",lty="dashed")
points(abscissa,truecurve$X0vera,type="l")
#pointwise interval: fix abscissa and see the interval for the relative fix abscissa

# Oversmoothing --> smoothing the data too much
# undersmoothing --> Use a basis system too rich
# generalized cross-validation (for balance undersmoothing and oversmoothing)
nbasis <- 6:30
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(0,1), nbasis[i], m)
  gcv[i] <- smooth.basis(abscissa, Xobs0, basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)] #min of gcv 

### SMOOTHING SPLINES --> reach basis system, retrieve smoothnessby penalization
breaks <- abscissa[((0:50)*2)+1]
basis <- create.bspline.basis(breaks, norder=m)
functionalPar <- fdPar(fdobj=basis, Lfdobj=3, lambda=1e-8)#basis, order of the derivative to be penalized, smoothing parameter
Xss <- smooth.basis(abscissa, Xobs0, functionalPar)
Xss0 <- eval.fd(abscissa, Xss$fd, Lfd=0)
Xss1 <- eval.fd(abscissa, Xss$fd, Lfd=1)
Xss2 <- eval.fd(abscissa, Xss$fd, Lfd=2)
df <- Xss$df   #  the degrees of freedom in the smoothing curve--> related to the penalization
#changin lambda, if df decrease --> stronger constrain on 3rd derivative
# Recommendation: when choosing the smoothing parameter, look at the
# derivatives vs central finite differences

gcv <- Xss$gcv #  the value of the gcv statistic

par(mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xss0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xss1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xss2 ,type="l",col="blue",lwd=2)
plot(basis)

# generalized cross-validation: choosing lambda
lambda <- c(1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12)
gcv <- numeric(length(lambda))
for (i in 1:length(lambda)){
  functionalPar <- fdPar(fdobj=basis, Lfdobj=3, lambda=lambda[i])  
  gcv[i] <- smooth.basis(abscissa, Xobs0, functionalPar)$gcv
}

#Contrained functions
### LOCAL POLYNOMIAL REGRESSION
library(KernSmooth)
m <- 5           # order of the polynomial
degree <- m-1    # degree of the polynomial
bw <- 0.05 # bandwidth: too large -->oersmoothing
#                       too small -->undersmoothing
Xsm0 <- locpoly(abscissa, Xobs0, degree=degree,
                bandwidth=bw, gridsize=length(abscissa), 
                range.x=range(abscissa))#posso mettere drv=grado derivata
Xsm0 <- Xsm0$y

# SMOOTHING for positive curves: help(smooth.pos)
#Exponential trasfromation
# y_j = exp(w(t_j)) + e_j 
# f(t) = exp(w(t)) 

# The function w(t) is unconstrained
# The function f(t) is monotone increasing
# w(t) is modeled via a basis expansion numerical methods are 
#used to compute the coefficients of the basis expansion

# SMOOTHING for monotone curves:  smooth.monotone
# f(t) = int_(t_0)^t exp(w(u)) du
# y_j = b_0 + b_1 * f(t_j) + e_j

# The function w(t) is unconstrained
# The function f(t) is monotone increasing
# b_1>0 for monotone increasing functions
# b_1<0 for monotone decreasing functions
# b_0 is the value of the function at t_0

# w(t) is modeled via a basis expansion 
# numerical methods are used to compute the coefficients of the
# basis expansion, as well as b_0, b_1

norder <- 6
nbasis <- nage - 2 + norder 
wbasis <- create.bspline.basis(rangeval = ageRng, nbasis = nbasis, 
                               norder = norder, breaks = age)
Lfdobj <- 3 #penalize 3th derivative   
lambda <- 10^(-0.5)  
cvecf <- matrix(0, nbasis, ncasef) 
Wfd0 <- fd(coef = cvecf, basisobj = wbasis)
growfdPar <- fdPar(fdobj = Wfd0, Lfdobj = Lfdobj, lambda = lambda)
growthMon <- smooth.monotone(argvals = age, y = hgtf, WfdParobj = growfdPar)



###FDA: FPCA ----
library(fda)
##Fourier basis
#if I see periodicity in the data --> Fourier basis
time <- 1:365 
#Metto il tempo sulle righe!!
# Choice 1: we set a high dimensional basis (interpolating)
  # Pros: no loss of information
  # Cons: possible overfitting 
data_W <- CanadianWeather$dailyAv[,,1]
basis.1 <- create.fourier.basis(rangeval=c(0,365),nbasis=365)
data_W.fd.1 <- Data2fd(y = data_W,argvals = time,basisobj = basis.1)
plot.fd(data_W.fd.1)

# Choice 2: reduced dimensionality (we set a low dimensional basis)
# Pros: the data are much smoother and the measurement error is filtered
# Cons: I could have lost important information
basis.2 <- create.fourier.basis(rangeval=c(0,365),nbasis=21)
data_W.fd.2 <- Data2fd(y = data_W,argvals = time,basisobj = basis.2)
plot.fd(data_W.fd.2) #smooth estimation, no noise

# Choice 3: compromise between 1 and 2
basis.3 <- create.fourier.basis(rangeval=c(0,365),nbasis=109)
data_W.fd.3 <- Data2fd(y = data_W,argvals = time,basisobj = basis.3)
plot.fd(data_W.fd.3)

#first 3 coefficen
data_W.fd.1$coefs[1:3,1]

# estimate of the mean and of the covariance kernel
library(fields)
#mean
plot.fd(data_W.fd.1)
lines(mean.fd(data_W.fd.1),lwd=3)
# covariance
eval.1 <- eval.fd(time,data_W.fd.1)
image.plot(time,time,(cov(t(eval.1))[1:365,])) 

##B-spline basis 
basis <- create.bspline.basis(rangeval=c(0,350),nbasis=21)
data_L.fd <- Data2fd(y = data_L,argvals = time,basisobj = basis)
plot.fd(data_L.fd, main="B-splines")
# Estimate of the mean and of the covariance kernel
layout(cbind(1,2))
plot.fd(data_L.fd,xaxs='i')
lines(mean.fd(data_L.fd),lwd=2)
eval <- eval.fd(time,data_L.fd)
image.plot(time, time, (cov(t(eval))[1:51,]))

### FPCA
# interpolated data (Choice 1)
plot.fd(data_W.fd.1,ylab='temperature')
pca_W.1 <- pca.fd(data_W.fd.1,nharm=5,centerfns=TRUE) #mharms: n* of pc to compute
#screeplot: [1: n^ of basis-1]
plot(pca_W.1$values[1:35],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_W.1$values)[1:35]/sum(pca_W.1$values),xlab='j',ylab='CPV',ylim=c(0.8,1))

# first two FPCs
par(mfrow=c(1,3))
plot(pca_W.1$harmonics[1,],col=1,ylab='FPC1',ylim=c(-0.1,0.08))
abline(h=0,lty=2)
plot(pca_W.1$harmonics[2,],col=2,ylab='FPC2',ylim=c(-0.1,0.08))

par(mfrow=c(1,2))
plot.pca.fd(pca_W.1, nx=100, pointplot=TRUE, harm=c(1,2), expand=0, cycle=FALSE)

#variance explained along the first 5 functional principal components
pca_W.1$varprop

# scatter plot of the scores
par(mfrow=c(1,2))
plot(pca_W.1$scores[,1],pca_W.1$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
points(pca_W.1$scores[35,1],pca_W.1$scores[35,2],col=2, lwd=4)

plot(pca_W.1$scores[,1],pca_W.1$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2",xlim=c(-400,250))
text(pca_W.1$scores[,1],pca_W.1$scores[,2],dimnames(data_W)[[2]], cex=1)
#see outliers
layout(1)
matplot(eval.1,type='l')
lines(eval.1[,35],lwd=4, col=2) 



###FDA: K-mean alignment ----
library(fdakma)
# Without alignment (just to perform clustering), let's try with 3 clusters:
set.seed(4)
fdakma_example_noalign_0der <- kma(
  x=x, y0=y0, n.clust = 3, 
  warping.method = 'NOalignment', #affine: combination of shift e dilation
  similarity.method = 'd0.pearson',# d1.pearson
  # similarity computed as the cosine between the original curves (correlation)
  center.method = 'k-means'
  #,seeds = c(1,11,21) # you can give a little help to the algorithm...
)
kma.show.results(fdakma_example_noalign_0der)

table(fdakma_example_noalign_0der$labels,
      fdakma_example_noalign_1der$labels, dnn=c("0der", "1der"))


# Warpings applied to the original data, colored according to  
# membership to the 3 clusters obtained via k-means, no alignment,
# d0.pearson similarity

plot(x, type="n", xlim=c(min(x),max(x)), ylim=c(min(x),max(x)+2), xlab="abscissa", ylab="warping")
title("Alignment affinities")
for(i in 1:30)(
  abline(a=fdakma_example$shift[i],b=fdakma_example$dilation[i], 
         col=fdakma_example_noalign_0der$labels[i])
)  

# How to choose the number of clusters and the warping method
kma.compare_example_2 <- kma.compare (
  x=x, y0=y0, y1=y1, n.clust = 1:3, 
  warping.method = c("NOaligment","shift", "dilation","affine"), 
  similarity.method = 'd1.pearson',
  center.method = 'k-means', 
  seeds = c(1,21,30),
  plot.graph=1)













