#Holy Bible

######################### Preliminar ###########################################
# Recall: qqplot to verify (qualitatively) the Gaussian assumption on the 
# distribution generating sample (m100 è una variabile!)
qqnorm(m100) # quantile-quantile plot
qqline(m100, col='red') # theoretical line
# Recall: Shapiro-Wilk test to verify (quantitatively) the Gaussian assumption on the
# distribution generating sample
shapiro.test(m100)

# Comandi base regression
regression <- lm(m200 ~ m100)
regression

summary(regression)

coef(regression)
vcov(regression)
residuals(regression)
fitted(regression)

points(m100, fitted(regression))
#plotta i punti della regressione sopra lo scatterplot

# Confidence and prediction intervals (command predict)
newdata <- data.frame(m100=c(10,11,12))
pred_nd <- predict(regression, newdata)

IC_nd <- predict(regression, newdata, interval = 'confidence', level = .99)
IP_nd <- predict(regression, newdata, interval = 'prediction', level = .99)

#fa tutti i grafici da se sgravatissimo
x11()
par (mfrow=c(2,2))
plot(regression)

aneurysm[,5] <- factor(aneurysm[,5]) #rende la colonna un factor

########################## Graphs ##############################################

# Scatter plot e istogrammi
plot(m100,m200)  #m100 e m200 sono colonne
hist(m100, prob=T)
hist(m200, prob=T)

plot(Data) #se dataset 2D con 2 covariate numeriche, puoi fare direttamente così
pairs(record) #record è un dataset piu di 2D
#penso che per i dataset numerici pairs e plot siano equivalenti

#Boxplot
boxplot(tourists, las=2, col='gold') #tourists è il dataset
boxplot(scale(x=tourists,center = T, scale=F), las=2, col='gold')
#scale scala i dati, scelgo dove centrarli,scale f vuol dire che li traslo e basta
#scale t invece divide anche per la standard deviation

plot(feed, weight) #se fai plot(numerica,categorica) ti fa i boxplot
#della numerica (asse y) per ogni categorica (asse x) 

# matplot
matplot(t(aneurysm.geometry),type='l')
matplot(t(aneurysm.geometry),type='l',col=color.position)
#Plot the columns of one matrix against the columns of another 
# (which often is just a vector treated as 1-column matrix).

# Pie chart (no ordering of levels) (solo per categoriche)
x11()
pie(table(district),col=rainbow(length(levels(district))))

# Barplot (levels are ordered) (solo per categoriche)
x11()
barplot(table(district)/length(district))
#istogramma delle frequenze realtive
plot(district)   # barplot of absolute frequences

################# PCA, verifica della teoria ###################################
#X è il dataset
M <- colMeans(X)   #mean di una matrice , quindi un vettore
S <- cov(X)   #covariance matrix
# we compute the eigenvectors and eigenvalues
eigen(S)
# Note. eigen(S)$vectors returns a matrix whose 
#     columns are the eigenvectors of S
# eigen(S)$values mi da gli autovalori
var.gen <- det(S)
var.tot <- sum( diag(S) )
#varianza totale e varianza generalizzata

library(car)
ellipse(M, S, 1, add=T,lwd=3, col='red')

#grafici
abline(a = M[2] - eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]*M[1], b = eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1], lty = 2, col = 'dark red', lwd = 2)
abline(a = M[2] - eigen(S)$vectors[2,2]/eigen(S)$vectors[1,2]*M[1], b = eigen(S)$vectors[2,2]/eigen(S)$vectors[1,2], lty = 2, col = 'red', lwd = 2)

############################# PCA vera #########################################
pc.tourists <- princomp(tourists, scores=T)
summary(pc.tourists)

# To obtain the rows of the summary:
# standard deviation of the components
pc.tourists$sd
# proportion of variance explained by each PC
pc.tourists$sd^2/sum(pc.tourists$sd^2)
# cumulative proportion of explained variance
cumsum(pc.tourists$sd^2)/sum(pc.tourists$sd^2)  #ovviamente l'ultimo fara sempre 1

#explained variance

load.tour <- pc.tourists$loadings #estrae i loading:i loading sono i coefficienti applicati 
#alle variabili originarie per determinare le componenti principali.

# graphical representation of the loadings of the first 8 principal components
x11()
par(mfcol = c(4,2))
for(i in 1:8) barplot(load.tour[,i], ylim = c(-1, 1), main=paste("PC",i))

# We compute the standardized variables
tourists.sd <- scale(tourists)
tourists.sd <- data.frame(tourists.sd) #creo nuovo dataset standardizzato

# scores
scores.tourists <- pc.tourists$scores
x11()
plot(scores.tourists[,1:2])
abline(h=0, v=0, lty=2, col='grey')
#come i dati si distribuiscono rispetto alle 2 principal component

biplot(pc.tourists) # instead of plotting points, plotto il numero della riga del dato

######################## Gaussain analysis in multivariate data ################
#prima studio le singole componenti, come nel caso 1D: se non sono normali x1 e x2,
#di sicuro non lo sara il vettore x=(x1,x2)

hist(X[,1], prob=T, ylab='density', xlab='X.1', main='Histogram of X.1',ylim=c(0,0.45))
lines((-1000):1000 /100, dnorm((-1000):1000 /100,mean(X[,1]),sd(X[,1])), col='blue', lty=2)
#guardando gli istogrammi cerco di vedere se i dati hanno un andamento normale

qqnorm(X[,1], main='QQplot of X.1',xlab='theoretical quantiles', ylab='sample quantiles')
qqline(X[,1])
#poi uso il qqplot sulle singole componenti, come fossi nel caso 1D
#tutto cio ha poco senso perche sto checkando la normalita solo in una direzione qualunque

shapiro.test(X[,1])

# we look at the directions of the PCs
pc.X <- princomp(X,scores=T)
#faccio le stesse cose con pc.X$scores[,1] al posto di X[,1], sia grafici che shapiro

### Consider the squared Mahalanobis distances of the data from the (sample) mean
### and test if they are a sample from a chi-square distribution
# Theorem: if X~N(mu,Sigma) r.v. in R^p, det(Sigma)>0
#          then d2(X,mu)=(X-mu)'Sigma^-1(X-mu) ~ Chi-sq(p)

M <- colMeans(stiff)  #rinominate in seguito: x.mean 
S <- cov(stiff)       #x.cov
d2 <- matrix(mahalanobis(stiff, M, S))

windows()
par(mfrow=c(1,2))
hist(d2, prob=T)
lines(0:2000/100, dchisq(0:2000/100,4), col='blue', lty=2)
qqplot(qchisq(seq(0.5/30, 1 - 0.5/30 ,by = 1/30), df = 4), d2,  main='QQplot di d2')
abline(0, 1)

### we can perform a chi.sq goodness-of-fit test
d2.class <- cut(d2, qchisq((0:6)/6, df = 4))
d2.freq  <- table(d2.class)
chisq.test(x = d2.freq, p = rep(1/6, 6), simulate.p.value = T)

# test of all the directions simultaneously
mcshapiro.test(stiff)
#output:  Wmin   <- mshapiro.test(t(X))$stat   # min(W) for the given sample
#pvalue <- sum(W < Wmin)/sim  # proportion of min(W) more extreme than the observed Wmin
#devst  <- sqrt(pvalue*(1-pvalue)/sim) #deviazione standard
#sim bo dato in input, ignoralo
#se il p value è molto alto accetto la normalita, guarda solo lui
# Wmin lo vorrei vicino a 1, ma basta giuardare il pvalue

mcshapiro.test(cbind(BC.x,BC.y,BC.z,BC.w))
#funziona non solo con i dataset ma anche con i vettori di covariate

#################### se trovo un pvalue bassissimo? ############################
#prova a togliere gli outlier
stiff.noout <- stiff[which(d2<7.5),]
#creo un nuovo daatset dove cavo i punti con una d2 distance troppo grande

#Box cox trasformation
# For lambda<1: observations <1 are "spread", observations >1 are "shrinked"
# For lambda>1: observations <1 are "shrinked", observations >1 are "spread"
# We compute the optimal lambda of the univariate Box-Cox transformation
# Caso 1D
lambda.x <- powerTransform(x)  # (command powerTransform of library car)
# Transformed sample with the optimal lambda (command bcPower of library car)
bc.x <- bcPower(x, lambda.x$lambda) # it transforms the data of the first argument 
# through the Box-Cox transformation with lambda given as second argument

# Caso multivariato
lambda <- powerTransform(cbind(x,y))  
BC.x <- bcPower(x, lambda$lambda[1])
BC.y <- bcPower(y, lambda$lambda[2])

#plot della box cox trasformation
xx <- seq(0,7,by=0.01)

plot(xx,box_cox(x=xx,lambda=lambda$lambda[1]),col='red',lty=1,type='l',xlab=expression(x),ylab=expression(x[lambda]),ylim=c(-5,10),asp=1)
title(main='Box-Cox transformation')
lines(xx,box_cox(x=xx,lambda=lambda$lambda[2]),col='blue')
points(1,0,pch=19,col='black')
abline(a=-1,b=1,lty=2,col='grey')
legend('bottomright',c(expression(lambda[x]),expression(lambda[y]),expression(paste(lambda,'=1'))),col=c('red','blue','grey'),lty=c(1,1,1))

# nel caso high dimensionality box cox fa schifo: buoni i p value one at time dello 
# shapiro ma terrificante l mc shapiro

################# Tests and confidence regions #################################
# x è il mio dataset
x.mean   <- sapply(x,mean)
x.cov    <- cov(x)
x.invcov <- solve(x.cov) #sarebbe s alla meno 1, serve alle baraccate di teoria

# Test sulla media , con dati normali (mu0 e alpha dati)
n <- dim(x)[1]  #numero di righe
p <- dim(x)[2]  #numero di colonne
# T2 Statistics
x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
# Radius of the ellipsoid
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
# Test: 
x.T2 < cfr.fisher   # alpha e mu0 devono essere forniti
# Rejection region: {x.T2>cfr.fisher}  (we reject for large values of the T2 statistics)
# true accetto , false rifiuto
P <- 1-pf(x.T2*(n-p)/((n-1)*p), p, n-p)     # p-value

#se alpha = p value, confidence inetrval interseca la mu0, quindi se faccio ellipse
# con alpha = p-value prendo sul bordo mu0, per definizione di p value

# Region of rejection (centered in mu0)
plot(x, asp = 1)
ellipse(mu0, shape=x.cov/n, sqrt(cfr.fisher), col = 'blue', lty = 2, center.pch = 16)
# We add a red point in correspondence of the sample mean
points(x.mean[1], x.mean[2], pch = 16, col ='red', cex = 1.5)
# Confidence region (centered in x.mean)
ellipse(x.mean, x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 2, lwd=2, center.cex=1)
# comando in piu per aggiungere il rettangolo
rect(T2[1,1],T2[2,1],T2[1,3],T2[2,3], border='red', lwd=2)
#calcolo del confidence interval, sempre piu grande dell'ellisse
T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2

# Asymptotic test, cambiano le distribuzioni (non serve check sulla normalita)
x.T2A   <- n * (x.mean-mu0) %*%  x.invcov  %*% (x.mean-mu0)  #cambia la distribuzione 
cfr.chisq <- qchisq(1-alpha,p)
x.T2A < cfr.chisq   # true accetto , false rifiuto
PA <- 1-pchisq(x.T2A, p)   #pvalue

# le asymptotic confidence e rejection regions sono piu larghe di quelle ottenute
# coi dati normali

# Bonferroni
k <- p # number of intervals I want to compute (set in advance)
cfr.t <- qt(1-alpha/(2*k),n-1)
Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf   #rettangolino piu piccino
#plot, da aggiungere al grafico di prima
rect(Bf[1,1],Bf[2,1],Bf[1,3],Bf[2,3], border='orange', lwd=2)

# il rettangolo di bonferroni è più piccolo per dimensioni basse (tipo 2 covariate)
# all'aumentare delle dimensioni diventa piu grande degli IC

#Plot of the confidence region in more dimensions
#prima calcola T2, come gia visto
matplot(1:4,1:4,pch='',ylim=range(stiff),xlab='Variables',ylab='T2 for a component', 
        main='Simultaneous T2 conf. int. for the components')
for(i in 1:4)segments(i,T2[i,1],i,T2[i,3],lwd=3,col=i)
points(1:4, T2[,2], pch=16, col=1:4)
#ho plottato i simultaneous confidence interval
points(1:4, mu0, lwd=3, col='orange')  #aggiungo mu0 per vedere se è dentro

#stessa cosa con bonferroni
matplot(1:4,1:4,pch='',ylim=range(stiff),xlab='Variables',ylab='Bf for a component',
        main='Simultaneous Bf conf. int. for the components')
for(i in 1:4) segments(i,Bf[i,1],i,Bf[i,3],lwd=2,col=i)
points(1:4, Bf[,2], pch=16, col=1:4)
points(1:4, mu0, lwd=3, col='orange')

# Test sulla media di paired multivariate gaussian
# esempio: 2 lab fanno la stessa roba, i dati agree?

# we compute the sample of differences
D <- data.frame(DBOD=effluent[,1]-effluent[,3], DSS=effluent[,2]-effluent[,4]) 
D   # nuovo dataset con solo 2 covariate: diff tra lab1 e lab2

# comando carino nel caso 2D: D è il dataset, plotta lui
plot(D, asp=1, pch=19, main='Dataset of Differences')
abline(h=0, v=0, col='grey35')
points(0,0, pch=19, col='grey35')

# dopo lo shapiro, se i dati sono normali fai il test con mu0=(0,0), easy

# sometimes a data is inside the rectangle but outside the ellipse: 
# i have to find the worst direction, that is the one of the Hotelling T2 statistics

# => From the theory
#    - the maximum is realized (Hotelling T2-statistics)
D.T2
#    - the distribution of the maximum is known
#    - the direction along which the maximum is realized is known
worst <- D.invcov %*% (D.mean-delta.0)
worst <- worst/sqrt(sum(worst^2))
worst
# Angle with the x-axis:
theta.worst <- atan(worst[2]/worst[1])+pi
theta.worst
# Confidence interval along the worst direction:
IC.worst  <- c( D.mean %*% worst - sqrt(cfr.fisher*(t(worst)%*%D.cov%*%worst)/n),
                D.mean %*% worst,
                D.mean %*% worst + sqrt(cfr.fisher*(t(worst)%*%D.cov%*%worst)/n) )
IC.worst
delta.0%*%worst
(IC.worst[1] < delta.0%*%worst) & (delta.0%*%worst < IC.worst[2])   
# se esce falso rifiuta, h0 è fuori dal rettangolo della worst direction

# Extremes of IC.worst in the coordinate system (x,y): plot
x.min <- IC.worst[1]*worst
x.max <- IC.worst[3]*worst
m1.ort <- -worst[1]/worst[2]
q.min.ort <- x.min[2] - m1.ort*x.min[1]
q.max.ort <- x.max[2] - m1.ort*x.max[1]
abline(q.min.ort, m1.ort, col='forestgreen', lty=2,lwd=1)
abline(q.max.ort, m1.ort, col='forestgreen', lty=2,lwd=1)

m1=worst[2]/worst[1] # worst direction
abline(0, m1, col='grey35')
segments(x.min[1],x.min[2],x.max[1],x.max[2],lty=1,lwd=2, col='forestgreen')

# nel lab 6 fa tante robe per dimostrare che è vero

# Test for repeated measure (richiede normalita)
# we build one of the possible contrast matrices 
C <- matrix(c(-1, 1, 0, 0,
              -1, 0, 1, 0,
              -1, 0, 0, 1), 3, 4, byrow=T)
C
# here we are looking at the effects on the pressure
# after 8, 16 and 24 hours from the instant the drug was given

Md <- C %*% M   #M è la media
Sd <- C %*% S %*% t(C)  #S è la covarianza
Sdinv <- solve(Sd)

# fai il solito test con queste 3 matrici! puoi fare tutto usando queste 3 
# invece di x.mean x.cov x.invcov

# Test for the mean of two independent Gaussian populations
# t1 e t2 data frame indipendenti
n1 <- dim(t1)[1] 
n2 <- dim(t2)[1] 
p  <- dim(t1)[2] 
#dopo aver calcolato media e covar per each dataset
Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2)
# we compare the matrices
list(S1=t1.cov, S2=t2.cov, Spooled=Sp)
Spinv   <- solve(Sp)
#nuovi test:
T2 <- n1*n2/(n1+n2) * (t1.mean-t2.mean-delta.0) %*% Spinv %*% (t1.mean-t2.mean-delta.0)
cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
T2 < cfr.fisher # if TRUE: no statistical evidence to reject H0 at level alpha
P <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
P  

# RIVEDI POCO CHIARO
# Caso alternativo. 2 dataset diversi, "almond" e "non almond" 
p.hat <- (x.mean1*n1+x.mean2*n2)/(n1+n2) # extimate of the proportion: weightetd mean
# è la p della bernouli
x.var <- (p.hat*(1-p.hat))  #sample variance: bernoulli

#usa i vettori invece del ciclo for
z.i <- (x.mean1-x.mean2)/sqrt(x.var*(1/n1+1/n2))
p.i <- ifelse(z.i<0, 2*pnorm(z.i),2*(1-pnorm(z.i))) 

p.Bf <- p.adjust(p.i, method='bonferroni') #aggiusto il pvalue con bonferroni
p.BH <- p.adjust(p.i, method='BH') #aggiusto con il false discovery rate

######################### Anova and Manova ######################################
##### One-way ANOVA
### Model: weigth.ij = mu + tau.i + eps.ij; eps.ij~N(0,sigma^2)
### Test:
### H0: tau.1 = tau.2 = tau.3 = tau.4 = tau.5 = tau.6 = 0
# il trattamento non ha effetto: H0

barplot(tapply(weight, feed, mean)) #istogramma per vedere gia graficamente se
#le medie sono diverse, weight e feed sono le 2 covariate, una num una cat

# One way anova: una num osservata vs una cat
n       <- length(feed)      # total number of obs.
ng      <- table(feed)       # number of obs. in each group
treat   <- levels(feed)      # levels of the treatment (sarebbero i nomi)
g       <- length(treat)     # number of levels (i.e., of groups)

### verify the assumptions:
# 1) normality (univariate) in each group (6 tests)
Ps <- c(shapiro.test(weight[ feed==treat[1] ])$p,
        shapiro.test(weight[ feed==treat[2] ])$p,
        shapiro.test(weight[ feed==treat[3] ])$p,
        shapiro.test(weight[ feed==treat[4] ])$p,
        shapiro.test(weight[ feed==treat[5] ])$p,
        shapiro.test(weight[ feed==treat[6] ])$p) 
Ps  # p value alto => è normale
# 2) same covariance structure (= same sigma^2)
Var <- c(var(weight[ feed==treat[1] ]),
         var(weight[ feed==treat[2] ]),
         var(weight[ feed==treat[3] ]),
         var(weight[ feed==treat[4] ]),
         var(weight[ feed==treat[5] ]),
         var(weight[ feed==treat[6] ])) 
Var  #faccio la varianza di ogni gruppo separatamente, solo qualitativo, inutile

# test of homogeneity of variances (questo è il test vero)
# H0: sigma.1 = sigma.2 = sigma.3 = sigma.4 = sigma.5 = sigma.6 
# H1: there exist i,j s.t. sigma.i!=sigma.j
bartlett.test(weight, feed) #p value alto => stessa varianza in ogni gruppo

fit <- aov(weight ~ feed)  #comando per lanova
summary(fit)
### How to read the summary:
#              Df   Sum Sq      Mean Sq      F value     Pr(>F)    
#  treat      (g-1) SStreat  SStreat/(g-1)  Fstatistic  p-value [H0: tau.i=0 for every i]
#  Residuals  (n-g) SSres     SSres/(n-g)  
# p value basso => signifcativo, dunque la media varia tra i gruppi, il treat ha effetto

#quale level è il responsabile? n-1 comparisons con bonferroni
k <- g*(g-1)/2
alpha= 0.05

Mediag  <- tapply(weight, feed, mean)
SSres <- sum(residuals(fit)^2)
S <- SSres/(n-g)

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

x11(width = 14, height = 7)
par(mfrow=c(1,2))
plot(feed, weight, xlab='treat', ylab='weight', col = rainbow(6), las=2)
# la linea grigia è lo 0, cioe nessuna diff: tutti gli ic che 
# non includono lo 0, hanno evidenza di avere una media diversa
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

##### One-way MANOVA
# Data exploration
colore <- rep(rainbow(3), each = 50)  #utile a salvare i colori
#ogni gruppo deve avere sempre lo stesso colore: n=50, 3 numero di levels
#penso funzioni solo per le osservazioni sortate

i1 <- which(species.name=='setosa')
i2 <- which(species.name=='versicolor')
i3 <- which(species.name=='virginica')
# in i1,i2,i3 salvo gli indici dei dati nel gruppo 1,2,3
# cosi poi li estraggo direttamente 

# in ogni boxplot, un gruppo / level diverso, in funzione delle covariate
# un colore per ogni covariata, per poterle confrontare tra i gruppi cioe tra i vari boxplot
x11(width=13)
par(mfrow=c(1,3))
boxplot(iris4[i1,], main='SETOSA',     ylim=c(0,8), col = rainbow(4))
boxplot(iris4[i2,], main='VERSICOLOR', ylim=c(0,8), col = rainbow(4))
boxplot(iris4[i3,], main='VIRGINICA',  ylim=c(0,8), col = rainbow(4))
# cazzata usare gli stessi colori per le covariate

# in ogni boxplot, una covariata diversa, in funzione dei gruppi
# un colore per ogni gruppo, confronto le differenze dentro lo stesso boxplot
x11(width=13)
par(mfrow=c(1,4))
boxplot(iris4[,1]~species.name, main='Sepal Length', ylim=c(0,8), col = rainbow(3))
boxplot(iris4[,2]~species.name, main='Sepal Width', ylim=c(0,8), col = rainbow(3))
boxplot(iris4[,3]~species.name, main='Petal Length', ylim=c(0,8), col = rainbow(3))
boxplot(iris4[,4]~species.name, main='Petal Width', ylim=c(0,8), col = rainbow(3))

# Test serio
### Model: X.ij = mu + tau.i + eps.ij; eps.ij~N_p(0,Sigma), X.ij, mu, tau.i in R^4
### Test:
### H0: tau.1 = tau.2 = tau.3  = (0,0,0,0)'
### H0: The membership to an iris species hasn't any significant effect on the mean
###     of X.ij (in any direction of R^4) 
### H1: There exists at least one direction in R^4 along which at least two species
###     have some feature significantly different

n1 <- length(i1)
n2 <- length(i2)
n3 <- length(i3)
n  <- n1+n2+n3

g  <- length(levels(species.name))  #numbers of groups
p  <- 4    # numero di covariate

### Verify the assumptions:
# 1)  normality (multivariate) in each group (3 tests)
Ps <- NULL
for(i in 1:g)
  Ps <- c(Ps, mcshapiro.test(iris[get(paste('i',i, sep='')),1:4])$p) 
Ps    #se i p value so alti è normale

# 2) same covariance structure (= same covariance matrix Sigma)
S  <-  cov(iris4)
S1 <-  cov(iris4[i1,])
S2 <-  cov(iris4[i2,])
S3 <-  cov(iris4[i3,])
#no statistical test, se fa a occhio
round(S1,digits=1)
round(S2,digits=1)
round(S3,digits=1)
#visual plot di quello che ho fatto sopra
x11(width=21)
par(mfrow=c(1,3))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S3, col=heat.colors(100),main='Cov. S3', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
# non sembrano uguali, ma vado avanti lo stesso (IN ESAME DILLO CHE VEDI CHE NON 
# SONO UGUALI) per risolvere la cosa prova una trasformazione

# Note: We can verify the assumptions a posteriori on the residuals of 
#       the estimated model, ecco perche vado avanti

fit <- manova(as.matrix(iris4) ~ species.name)
summary.manova(fit,test="Wilks")  #specifica quale test statistic vuoi usare
# Exact tests for p<=2 or g<=3 already implemented in R: pvalue baaso => H1

summary.aov(fit)
# fa 4 test separati, una per covariata: per ogni covariata mi dice se 
# l'appartenenza a un gruppo è significativa per rifiutare, ma non so quali gruppi
# influenzano di più e rispetto a quale covariata => CI

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
# come al solito guardo quali intervalli contengono lo 0

x11()
par(mfrow=c(1,4))
mg <- rbind(m1,m2,m3)
sp.name <- c('Sepal Length','Sepal Width', 'Petal Length', 'Petal Width')
for(k in 1:4){
  plot(c(1,g*(g-1)/2),ylim=c(-4,4), xlim=c(1,3), pch='', 
       xlab='pairs treat', ylab=paste('CI tau',k), 
       main=paste('CI tau',sp.name[k]))
  lines (c(1,1), c(CI[[1]][k,1],CI[[1]][k,2])); 
  points(1, mg[1,k]-mg[2,k], pch=16); 
  points(1, CI[[1]][k,1], col=rainbow(g)[2], pch=16); 
  points(1, CI[[1]][k,2], col=rainbow(g)[1], pch=16);  
  lines (c(2,2), c(CI[[2]][k,1],CI[[2]][k,2])); 
  points(2, mg[1,k]-mg[3,k], pch=16);
  points(2, CI[[2]][k,1], col=rainbow(g)[3], pch=16); 
  points(2, CI[[2]][k,2], col=rainbow(g)[1], pch=16);
  lines (c(3,3), c(CI[[3]][k,1],CI[[3]][k,2])); 
  points(3, mg[2,k]-mg[3,k], pch=16);
  points(3, CI[[3]][k,1], col=rainbow(g)[3], pch=16); 
  points(3, CI[[3]][k,2], col=rainbow(g)[2], pch=16);  
  abline(h=0)
}

##### Two-ways ANOVA
##### (p=1, g=2, b=2) g: #levels, p: #covariate, b: #ways, cioe 2 factor

# devi saper creare i dati in esame
# NON HA SENSO CONTROLLARE LA NORMALITA SU COSI POCHI DATI, mi tocca assumerla
km          <- c(18.7, 16.8, 20.1, 22.4, 14.0, 15.2, 22.0, 23.3)
distr       <- factor(c('Esso','Esso','Esso','Esso','Shell','Shell','Shell','Shell'))
benz        <- factor(c('95','95','98','98','95','95','98','98'))
distr_benz  <- factor(c('Esso95','Esso95','Esso98','Esso98','Shell95','Shell95','Shell98','Shell98'))
#distr benz laggiungo io: sarebbe la categorica "prodotto" 
g <- length(levels(distr))
b <- length(levels(benz))
n <- length(km)/(g*b)

M           <- mean(km)
Mdistr      <- tapply(km,      distr, mean)
Mbenz       <- tapply(km,       benz, mean)
Mdistr_benz <- tapply(km, distr_benz, mean)

### Model with interaction (complete model): 
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N(0,sigma^2), 
###     i=1,2 (effect station), j=1,2 (effect gasoline)

fit.aov2.int <- aov(km ~ distr + benz + distr:benz) #come fare aov su vettori 
# e non su un dataset, uso i : per usare 2 cate, conta anche l'interazione tra loro
# in pratica anova con 3 categoriche diverse, le 2 che gli do e la loro inetraction
# infatti nel summary ho 3 test
fit.aov2.int <- aov(km ~ distr*benz)  #notazione equivalente
summary.aov(fit.aov2.int)

### Test:    
### 1) H0: gamma.11 = gamma.12 = gamma.21 = gamma.22 = 0    vs   H1: (H0)^c 
###    H0: There is no significant interaction between the factors station
###        and gasoline in terms of performances
###    H1: There exists a significant interaction between the factors station 
###        and gasoline in terms of performances
###
### 2) H0: tau.1 = tau.2 = 0    vs   H1: (H0)^c   
###    H0: The effect "gas station" doesn't significantly influence performances 
###    H1: The effect "gas station" significantly influences performances
###
### 3) H0: beta.1 = beta.2 = 0    vs   H1: (H0)^c
###    H0: The effect "gasoline" doesn't significantly influence performances
###    H1: The effect "gasoline" significantly influences performances

### Additive model: non conto l'interazione se non è molto significativa
### X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N(0,sigma^2), 
###     i=1,2 (effect station), j=1,2 (effect gasoline)
fit.aov2.ad <- aov(km ~ distr + benz)
summary.aov(fit.aov2.ad)
# Remark: by removing the interaction, the residual degrees of freedom increase! 
# qui ho solo test 2 e test 3

### Example: global test for the significance of the two treatments 
###          (model without interaction) (Test 2)
SSdistr <- sum(n*b*(Mdistr - M)^2)              # or from the summary: 1.53    
SSbenz  <- sum(n*g*(Mbenz  - M)^2)              # or from the summary: 66.70
SSres   <- sum((km - M)^2) - (SSdistr+SSbenz)   # or from the summary: 16.37

Ftot      <- ( (SSdistr + SSbenz) / ((g-1)+(b-1)))/(SSres / (n*g*b-g-b+1))
Ptot      <- 1 - pf(Ftot, (g-1)+(b-1), n*g*b-g-b+1) # attention to the dof!
Ptot   # il treatment è significativo, ma non so per colpa di chi, porcata

### Reduced additive model (ANOVA one-way, b=2): 
### X.jk = mu + beta.j + eps.jk; eps.jk~N(0,sigma^2), j=1,2 (effect gasoline)
fit.aov1 <- aov(km ~ benz)
summary.aov(fit.aov1)

SSres <- sum(residuals(fit.aov1)^2)
### Interval at 90% for the differences (reduced additive model)
### [b=2, thus one interval only]
IC <- c(diff(Mbenz) - qt(0.95, (n*g-1)*b) * sqrt(SSres/((n*g-1)*b) *(1/(n*g) + 1/(n*g))), 
        diff(Mbenz) + qt(0.95, (n*g-1)*b) * sqrt(SSres/((n*g-1)*b) *(1/(n*g) + 1/(n*g))))
names(IC) <- c('Inf', 'Sup')
IC    # IC for mu(98)-mu(95), se non ce lo 0, differenza tra le 2 medie consisntente

# SE FAI ONE WAY ANOVA VERIFICA LE IPOTESI

##### Two-ways MANOVA, (p=3, g=2, b=2)
# caso con un dataset 
plastic <- read.table('T6-4.dat',col.names=c('Ex','Ad','Tr','Gl','Op'))
Ex   <- factor(plastic$Ex, labels=c('L','H')) # Treat.1
Ad   <- factor(plastic$Ad, labels=c('L','H')) # Treat.2
#costruisco la covariata prodotto
ExAd <- Ex
levels(ExAd) <- c('LL','LH','HL','HH')
ExAd[Ex=='L' & Ad=='L'] <- 'LL'
ExAd[Ex=='L' & Ad=='H'] <- 'LH'
ExAd[Ex=='H' & Ad=='L'] <- 'HL'
ExAd[Ex=='H' & Ad=='H'] <- 'HH'

#Data exploration per plot cute
boxplot(plastic3[,1]~ExAd) 
#boxplot(plastic3[,feature numerica]~feature categotica

#test normalita e hom della cov => vedi anova (scorri i levels usando Exad)

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

# ANOVA on the components (we look at the 3 axes-directions in R^3 separately)
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

qT <- qt(1 - alpha / (2 * k), g*b*n-g-b+1)
# the degrees of freedom of the residuals on the additive model are
# g*b*n-g-b+1

mExL  <- sapply(plastic3[Ex=='L',],mean)
mExH  <- sapply(plastic3[Ex=='H',],mean)
infEx <- mExH-mExL - qT * sqrt( diag(W)/(g*b*n-g-b+1) * (1/10+1/10) ) #estremo inferiore
supEx <- mExH-mExL + qT * sqrt( diag(W)/(g*b*n-g-b+1) * (1/10+1/10) ) #estremo superiore

mAdL  <- sapply(plastic3[Ad=='L',],mean)
mAdH  <- sapply(plastic3[Ad=='H',],mean)
infAd <- mAdH-mAdL - qT * sqrt( diag(W)/(g*b*n-g-b+1) * (1/10+1/10) )
supAd <- mAdH-mAdL + qT * sqrt( diag(W)/(g*b*n-g-b+1) * (1/10+1/10) )

IC2   <- list(ExH_ExL=cbind(infEx, supEx), AdH_AdL=cbind(infAd, supAd))
IC2

#################### Classifiers #############################################
### Linear and Quadratic Discriminant Analysis
### Example 1 (2 classes, univariate), LDA (univariate)

A <- which(group=='A')   # Group A
B <- which(group=='B')   # Group B

# Assumptions:
# 1) if L=i, X.i ~ N(mu.i, sigma.i^2), i=A,B
# 2) sigma.A=sigma.B
# 3) c(A|B)=c(B|A) (equal misclassification costs)

# 1) normality (univariate) within the groups
shapiro.test(cyto[A,1])
shapiro.test(cyto[B,1])

# 2) equal variance (univariate)
var.test(cyto[A,1],cyto[B,1]) # p value alto => accetto
#lo puoi usare solo se hai una sola covariata, sennò barlett

library(MASS)
help(lda)

cyto.lda <- lda(group ~ Infg)
cyto.lda

# posterior probability and classification for x=0
x <- data.frame(Infg = 0)  #creo un data frame con un solo valore, che sarebbe 0
# The command predict() returns a list containing (see the help of predict.lda):
# - the class associated with the highest posterior probability 
predict(cyto.lda, x)$class
# - the posterior probabilities for the classes
predict(cyto.lda, x)$posterior
# - in lda: the coordinates of the canonical analysis of Fisher
#           (Fisher's discriminant scores)
predict(cyto.lda, x)$x

#graph
x11()
plot(Infg[A], rep(0, length(A)), pch=16, col='blue', ylim=c(0,1),
     xlab='x', ylab='estimated posterior', main="LDA", xlim=range(Infg))
points(Infg[B], rep(0, length(B)), pch=16, col='red')
abline(v=0, col='grey')

# posterior probability for a grid of x's
x <- data.frame(Infg=seq(-10, 35, 0.5))   # come data frame ora ho un vettore
cyto.LDA.A <- predict(cyto.lda, x)$posterior[,1] # posterior probability for class A
cyto.LDA.B <- predict(cyto.lda, x)$posterior[,2] # posterior probability for class B
#ho salvato la colonna 1, cioe le prob di stare in A, same for B
predict(cyto.lda, x)$class #le classi predictate

lines(x[,1], cyto.LDA.A, type='l', col='blue', xlab='x', ylab='estimated posterior', main="LDA")
lines(x[,1], cyto.LDA.B, type='l', col='red')
abline(h = 0.5)
legend(-10, 0.9, legend=c('P(A|X=x)', 'P(B|X=x)'), fill=c('blue','red'), cex = 0.7)

# set prior probabilities
cyto.lda.1 <- lda(group ~ Infg, prior=c(0.95,0.05))  #basta aggiungere prior ecc...

### k-nearest neighbor classifier
library(class)
help(knn)  #usalo per vedere cosa è ogni parametro che gli do

cyto.knn <- knn(train = Infg, test = x, cl = group, k = 3, prob=T)
cyto.knn.class <- (cyto.knn == 'B')+0    #trasforma le var in 0 e 1
cyto.knn.B <- ifelse(cyto.knn.class==1, 
                     attributes(cyto.knn)$prob, 
                     1 - attributes(cyto.knn)$prob) #salvo le prob, in particolare:
#save the vote of the knn alg: cosa vota ogni nodo, se vedo 1 tutti e 3 votano per B
#0 tutti per A ecc, se supero 0.5 predicto B => ottengo una step function 

x11()
plot(x[,1], cyto.LDA.B, type='l', col='red', lty=2, xlab='x', ylab='estimated posterior')
points(x[,1], cyto.knn.B, type='l', col='black', lty=1)
legend(-10, 0.75, legend=c('LDA','knn'), lty=c(2,1), col=c('red','black'))
#plot rapidi di lda e knn, per la knn prova piu k

### Example 2 (3 classes, bivariate)
#uso iris: i1,i2,i3,g,n1,n2,n3,p riga 479 in poi
# Jittering
set.seed(1)
iris2 <- iris2 + cbind(rnorm(150, sd=0.025))    # jittering, aggiungo rumore
# per evitare overlap nei plot, utile per vedere se la mia analisi è robust

# Linear Discriminant Analysis (LDA), caso multivariate
lda.iris <- lda(iris2, species.name)
lda.iris

# "coefficients of linear discriminants" and "proportion of trace":
# Fisher discriminant analysis. 
# In particular:
# - coefficients of linear discriminants: versors of the canonical directions
#   [to be read column-wise]
# - proportion of trace: proportion of variance explained by the corresponding 
#   canonical direction

Lda.iris <- predict(lda.iris, iris2) #modo per evaluate il classifiers

# Estimate of AER (actual error rate):
# 1) APER (apparent error rate)
# 2) estimate of AER by cross-validation

# 1) Compute the APER
#Lda.iris$class : assigned classes , species.name : true labels
table(class.true=species.name, class.assigned=Lda.iris$class)
errors <- (Lda.iris$class != species.name)
APER   <- sum(errors)/length(species.name) 
#valori fuori dalla diagonale in table/totale, molto ottimistico, testo sul training set

# 2) Compute the estimate of the AER by leave-one-out cross-validation
LdaCV.iris <- lda(iris2, species.name, CV=TRUE)  # specify the argument CV
table(class.true=species.name, class.assignedCV=LdaCV.iris$class)
errorsCV <- (LdaCV.iris$class != species.name)
AERCV   <- sum(errorsCV)/length(species.name)
# Remark: correct only if we estimate the priors through the sample frequencies!

# Plot the partition induced by LDA

x11()
plot(iris2, main='Iris Sepal', xlab='Sepal.Length', ylab='Sepal.Width', pch=20)
points(iris2[i1,], col='red', pch=20)
points(iris2[i2,], col='green', pch=20)
points(iris2[i3,], col='blue', pch=20)
legend("topright", legend=levels(species.name), fill=c('red','green','blue'), cex=.7)
points(lda.iris$means, pch=4,col=c('red','green','blue') , lwd=2, cex=1.5)
#plotto con le x le medie

x  <- seq(min(iris[,1]), max(iris[,1]), length=200) #sorto la prima colonna
y  <- seq(min(iris[,2]), max(iris[,2]), length=200)
xy <- expand.grid(Sepal.Length=x, Sepal.Width=y)

z  <- predict(lda.iris, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - pmax(z[,2], z[,3])  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- z[,2] - pmax(z[,1], z[,3])  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    
z3 <- z[,3] - pmax(z[,1], z[,2])  # P_3*f_3(x,y)-max{P_j*f_j(x,y)}

# Plot the contour line of level (levels=0) of z1, z2, z3: 
# P_i*f_i(x,y)-max{P_j*f_j(x,y)}=0 i.e., boundary between R.i and R.j 
# where j realizes the max.
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)
contour(x, y, matrix(z3, 200), levels=0, drawlabels=F, add=T)

### Quadratic Discriminand Analysis (QDA)
###---------------------------------------
#qda invece che lda, la funzione è la stessa: stessi comandi

### knn-classifier
###----------------
# Plot the partition induced by knn, mi da piu flessibilità ma meno interpretabilità
#x11() bla bla bla come sopra
x  <- seq(min(iris[,1]), max(iris[,1]), length=200)
y  <- seq(min(iris[,2]), max(iris[,2]), length=200)
xy <- expand.grid(Sepal.Length=x, Sepal.Width=y)

iris.knn <- knn(train = iris2, test = xy, cl = iris$Species, k = k)
z  <- as.numeric(iris.knn)
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5), drawlabels=F, add=T)

##### Problem 2 of 9/09/2009
# utility: come costruire unico dataset e variabile factor
true <- read.table('moneytrue.txt',header=TRUE)
false <- read.table('moneyfalse.txt',header=TRUE)
banknotes <- rbind(true,false)  #costruisco un unico dataset
vf <- factor(rep(c('true','false'),each=100), levels=c('true','false'))

# misclassification costs
c.vf <- 10
c.fv <- 0.05
#prior probabilities
pf <- 0.001
pt <- 1-0.001
# Prior modified to account for the misclassification costs
prior.c <- c(pt*c.fv/(c.vf*pf+c.fv*pt), pf*c.vf/(c.vf*pf+c.fv*pt))

###-------------------------------------------------------------------------
### FISHER DISCRIMINANT ANALYSIS
###-----------------------------------
### Let's change viewpoint: we look for the directions that highlight
### the discrimination among groups -> we look for the canonical directions
# riprendo iris
lda.iris <- lda(iris2, species.name) # me le dice lui
lda.iris$scaling #directions
a1=lda.iris$scaling[,1]
a2=lda.iris$scaling[,2]

### How are the data classified?
# Compute the canonical coordinates of the data
cc1.iris <- as.matrix(iris2)%*%a1
cc2.iris <- as.matrix(iris2)%*%a2
#salvo le coordinate
coord.cc <- cbind(cc1.iris,cc2.iris)
# Compute the coordinates of the mean within groups along the canonical directions
cc.m1 <- c(m1%*%a1, m1%*%a2)
cc.m2 <- c(m2%*%a1, m2%*%a2)
cc.m3 <- c(m3%*%a1, m3%*%a2)

# Assign data to groups
f.class=rep(0, n)
for(i in 1:n) # for each datum
{
  # Compute the Euclidean distance of the i-th datum from mean within the groups
  dist.m=c(d1=sqrt(sum((coord.cc[i,]-cc.m1)^2)),
           d2=sqrt(sum((coord.cc[i,]-cc.m2)^2)),
           d3=sqrt(sum((coord.cc[i,]-cc.m3)^2)))
  # Assign the datum to the group whose mean is the nearest
  f.class[i]=which.min(dist.m) # sarebbe il nuovo class assigned
}

### How do I classify a new observation?
x.new <- c(5.85, 2.90)
# compute the canonical coordinates
cc.new <- c(x.new%*%a1, x.new%*%a2)
# compute the distance from the means
dist.m <- c(d1=sqrt(sum((cc.new-cc.m1)^2)),
            d2=sqrt(sum((cc.new-cc.m2)^2)),
            d3=sqrt(sum((cc.new-cc.m3)^2)))
# assign to the nearest mean
which.min(dist.m) # stampa il gruppo

# riga 706 lab 8 da li in poi vari grafici

### Support Vector Machines
# cerco un sottopsazio che divida lo spazio 
# Fit the Support Vector Classifier (kernel = "linear") given a cost C
dat <- data.frame(x=x, y=as.factor (y)) # è il mio dataset
svmfit <- svm(y~., data=dat , kernel ='linear', cost =10, scale =FALSE )
summary(svmfit)

x11()
par(mfrow=c(1,2))
plot(svmfit , dat, col =c('salmon', 'light blue'), pch=19, asp=1)

# support vectors are indicated with crosses, i tondi sono i dati
# they are: (estrae i support vector)
svmfit$index

### altri 2 plot, piu carini linea 864 lab 8

# If we try to change the cost parameter we get more support points
# (higher bias, lower variance), se c piu basso ho piu vectors

# To set the parameter C we can use the function tune(),
# which is based on cross-validation (10-fold)
set.seed (1)
tune.out <- tune(svm,y~.,data=dat ,kernel = 'linear',
                 ranges =list(cost=c(0.001 , 0.01, 0.1, 1,5,10,100) ))
summary(tune.out)
# Extract the best model from the result of tune
bestmod <- tune.out$best.model
summary(bestmod)

# Prediction for a new observation (command predict())
# given testdat
ypred <- predict(bestmod,testdat)
table(true.label=testdat$y, assigned.label =ypred )

# If the classes are separable, setting a high value for the cost function
# leads to the maximal margin classifier (i.e., it returns the classification
# provided by the best separating hyperplane)

### Non-linear case, classi non separabili
# Randomly split in train and test
train <- sample (200 ,100)
# Fit a Support Vector Machine (kernel = "radial") given a cost C
svmfit <- svm(y~., data=dat [train ,], kernel ='radial', gamma =1, cost =1)
summary(svmfit)

# Misclassification error on the training set
table(true=dat[train ,"y"], pred=predict (svmfit ,
                                          newdata =dat[train ,]))
# Misclassification error on the test set
table(true=dat[-train ,"y"], pred=predict (svmfit ,
                                           newdata =dat[-train ,]))
# Increasing the cost decreases the errors on the training set,
# at the expense of a more irregular boundary e errore piu alto nel test set

# Set parameters via CV: sta volta setto pure gamma !
tune.out <- tune(svm , y~., data=dat[train ,], kernel ='radial',
                 ranges =list(cost=c(0.1 ,1 ,10 ,100 ,1000),
                              gamma=c(0.5,1,2,3,4) ))

###################### Unsupervised Learning ###########################
# studio iris senza le labels

# mischia sempre i dati all'inizio!
# actually, the data are never ordered according to (unknown) labels
misc <- sample(150)
iris4 <- iris4[misc,]

# compute the dissimilarity matrix of the data
# we choose the Euclidean metric (and then we look at other metrics)
iris.e <- dist(iris4, method='euclidean')
# with other metrics:
iris.m <- dist(iris4, method='manhattan')
iris.c <- dist(iris4, method='canberra')

x11()  #confronto grafico delle distanze
par(mfrow=c(1,3))
image(1:150,1:150,as.matrix(iris.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j' )
image(1:150,1:150,as.matrix(iris.c), main='metrics: Canberra', asp=1, xlab='i', ylab='j' )
image(1:150,1:150,as.matrix(iris.m), main='metrics: Manhattan', asp=1, xlab='i', ylab='j' )

# hierarical clustering
iris.es <- hclust(iris.e, method='single')
iris.ea <- hclust(iris.e, method='average')
iris.ec <- hclust(iris.e, method='complete')
# if we want more detailed information on euclidean-complete clustering:
names(iris.ec)
iris.ec$merge  # order of aggregation of statistical units / clusters
iris.ec$height # distance at which we have aggregations
iris.ec$order  # ordering that allows to avoid intersections in the dendrogram

# plot of the dendrograms
x11()
par(mfrow=c(1,3))
plot(iris.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(iris.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(iris.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.es, k=2) #plotta rettangoli rossi nel dendogram

# Fix k=2 clusters:
cluster.ec <- cutree(iris.ec, k=2) # euclidean-complete:

# interpret the clusters
table(label.true = species.name[misc], label.cluster = cluster.es)

coph.es <- cophenetic(iris.es) # compute the cophenetic matrices
es <- cor(iris.e, coph.es) # compute cophenetic coefficients, se alto è meglio

# compare with dissimilarity matrix (Euclidean distance)
x11()
par(mfrow=c(1,2))
image(as.matrix(iris.e), main='Euclidean', asp=1 )
image(as.matrix(coph.es), main='Single', asp=1 )

# single linkage fallisce con clusters circolari, se ho palle di dati meglio gli altri
# single linkage vince se hai un andamento lineare
# all'esame plottali tutti contro e confrontali
x11()
par(mfrow=c(1,3))  #X è il dataset
plot(X, xlab='Var 1', ylab='Var 2', main = 'Single linkage', col=ifelse(cluster.es==1,'red','blue'), pch=16, asp=1)
plot(X, xlab='Var 1', ylab='Var 2', main = 'Average linkage', col=ifelse(cluster.ea==1,'red','blue'), pch=16, asp=1)
plot(X, xlab='Var 1', ylab='Var 2', main = 'Complete linkage', col=ifelse(cluster.ec==1,'red','blue'), pch=16, asp=1)

clustw <- hclust(d, method='ward.D2') #ward linkage, altro tipo di clustering

##### K-means method
result.k <- kmeans(Q, centers=2) # Centers: fixed number of clusters,lo scegli all'inizio

result.k$cluster      # labels of clusters
result.k$centers      # centers of the clusters
result.k$totss        # tot. sum of squares
result.k$withinss     # sum of squares within clusters
result.k$tot.withinss # sum(sum of squares within cluster)
result.k$betweenss    # sum of squares between clusters
result.k$size         # dimension of the clusters

x11()
plot(Q, col = result.k$cluster+1)

open3d()
plot3d(Q, size=3, col=result.k$cluster+1, aspect = F) 
points3d(result.k$centers,size=10)

### How to choose k:
### 1) evaluate the variability between the groups with respect to 
###   the variability withing the groups 
### 2) evaluate the result of hierarchical clustering (not recommended,
###    quite computationally expensive)
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

############################## Linear models ###################################
### Model:
### distance = beta_0 + beta_1 * speed + beta_2 * speed^2 + Eps
### Assumptions:
## 1) Parameter estimation: E(Eps) = 0  and  Var(Eps) = sigma^2 
## 2) Inference:            Eps ~ N(0, sigma^2)
##### 1) Estimate of the parameters
##### Assumptions: E(Eps) = 0  and  Var(Eps) = sigma^2 
#####--------------------------------------------------
fm <- lm(distance ~ speed1 + speed2) # distance e speed sono colonne
summary(fm)  # in alternativa puoi usare cars$dist ecc..

fitted(fm)        # y hat
residuals(fm)     # eps hat
coefficients(fm)  # beta_i
vcov(fm)          # cov(beta_i) , covariance matrix dei coefficients
fm$rank # order of the model [r+1]
fm$df   # degrees of freedom of the residuals [n-(r+1)]
hatvalues(fm) # h_ii
rstandard(fm) # standardized residuals
sum(residuals(fm)^2)/fm$df  # estimate of sigma^2

##### Inference on the parameters
##### Assumption: Eps ~ N(0, sigma^2)
#####-------------------------------------------
### Test (Fisher):
# H0: (beta1, beta2) == (0, 0) vs H1: (beta1, beta2) != (0, 0)
linearHypothesis(fm, rbind(c(0,1,0), c(0,0,1)), c(0,0)) 
# i vettori sono i coeff che metto prima dei parametri: tipo sopra sto 
# confrontando b1 (primo vettore) con b2 (secondo)
summary(fm)
# i test li puoi fare a mano usando la t statistics (solo Oat!)
# linearHypothesis serve per testarli both contemporaneamente

p <- 2  # number of tested coefficients
r <- 2  # number of regressors

# Confidence region:
c(coefficients(fm)[2], coefficients(fm)[3])  # center: point estimate
eigen(vcov(fm)[2:3,2:3])$vectors  # Direction of the axes?

x11()
plot(coefficients(fm)[2], coefficients(fm)[3], xlim = c(-6,6), ylim = c(-6,6), asp=1, xlab='beta1', ylab='beta2')
ellipse(coefficients(fm)[2:3], vcov(fm)[2:3,2:3], sqrt(p*qf(1-0.05,p,n-(r+1))))
abline(v=0)
abline(h=0)

confint(fm, level= 1-0.05/p)[2:3,]  # Bonferroni correction!
# Note: confint() returns the confidence intervals one-at-a-time;
# to have a global level 95% we need to include a correction

### altro esempio di Test:
# H0: (beta0+beta2, beta1) == (0,0) vs H1: (beta0+beta2, beta1) != (0,0)
C <- rbind(c(1,0,1), c(0,1,0))
linearHypothesis(fm, C, c(0,0))

##### Confidence intervals for the mean & prediction (new obs)
##### Assumption: Eps ~ N(0, sigma^2)
#####---------------------------------
# Command predict()
Z0.new <- data.frame(speed1=10, speed2=10^2) # questo è il nuovo dato
# Conf. int. for the mean
Conf <- predict(fm, Z0.new, interval='confidence', level=1-0.05)  
Conf  # intervallo sulla media campionaria delle nuove oobs
# Pred. int. for a new obs
Pred <- predict(fm, Z0.new, interval='prediction', level=1-0.05)  
Pred  # intervallo sulle nuove obs, piu grande

#stessa cosa, ma su un dataset z0
Z0   <- data.frame(cbind(speed1=seq(0, 30, length=100), 
                         speed2=seq(0, 30, length=100)^2))
Conf <- predict(fm, Z0, interval='confidence')
Pred <- predict(fm, Z0, interval='prediction')

x11()
plot(cars, xlab='Speed', ylab='Stopping distance', las=1, xlim=c(0,30), ylim=c(-45,150))
lines(Z0[,1], Conf[,'fit'])
lines(Z0[,1], Conf[,'lwr'], lty=2, col='red', lwd=2)
lines(Z0[,1], Conf[,'upr'], lty=2, col='red', lwd=2)
lines(Z0[,1], Pred[,'lwr'], lty=3, col='gold', lwd=2)
lines(Z0[,1], Pred[,'upr'], lty=3, col='gold', lwd=2)
# vicino al baricentro ho meno variability, previsioni migliori

#in alternativa usa matplot:
matplot(Z0,Conf,add=T,type='l',col=c('black','blue','blue'),lwd=2,lty=2)
matplot(Z0,Pred,add=T,type='l',col=c('black','green','green'),lwd=2,lty=2)
legend('topleft', legend=c('regression line','confidence intervals','prediction intervals'),
       col=c('black','blue','green'), lwd=2, cex=0.85)

##### Verify assumptions
#####--------------------
x11()
par(mfrow=c(2,2))
plot(fm)  # fa tutti i plot utili
shapiro.test(residuals(fm))
# If the pvalue of the Shapiro-Wilk Test is greater than 0.05, the data is normal

#plot della retta di regressione nel caso 2D
result <- lm(Y ~ X)
coef <- result$coef
plot(X, Y)
abline(coef[1],coef[2], lwd=2,col='red')

### Example 4: Earthquakes
# Model:
# depth = beta0 + beta1*lat + beta2*long + eps
fit  <- lm(depth ~ lat + long, data=Q) # modo per chiamare direttamente dal dataset Q

### Model 2 (dummy variable)
###--------------------------
dummy <- clusterw - 1   # 0 = red, clusterw arriva da cutree
# 1 = green
Qd <- cbind(Q, dummy)

# Model:
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
# test: are the two planes needed? (sto facendo il test senza gaussianity!)
A <- rbind(c(0,0,0,1,0,0), c(0,0,0,0,1,0), c(0,0,0,0,0,1))
b <- c(0,0,0)
linearHypothesis(fitd, A, b)
# sto testando se la categorical variable è necessaria. pvalue piccolo quindi si

# skippata riduzione modello, semplicemente togli lat:dummy
# come fare grafico 3D:
open3d()
par3d(windowRect=c(680,40,1350,720))
points3d(x=Q$lat, y=Q$long, z=Q$depth, size=4, col=clusterw+1, aspect = T)
axes3d()
points3d(x=Qd$lat, y=Qd$long, z=fitted(fitD), size=4, col = 'blue')
surface3d(range(Q$lat), range(Q$long), 
          matrix(predict(fitD, expand.grid(lat = range(Q$lat), long = range(Q$long), dummy=c(1,1))),2,2),
          alpha = 0.5, col='green')
surface3d(range(Q$lat), range(Q$long), 
          matrix(predict(fitD, expand.grid(lat = range(Q$lat), long = range(Q$long), dummy=c(0,0))),2,2),
          alpha = 0.5, col='red')
