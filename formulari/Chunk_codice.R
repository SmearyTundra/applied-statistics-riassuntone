# 1)COMANDI GEN ##################################################################

#tilde ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Per campioni gaussiani 1-dim:
# Xn ~ N(mu, sigma^2/n)
# S ~ sigma^2*chi-sq(n-1)/n-1

cat("\014") #to clear the console
dataset = read.table('nome.txt', header=T) 
dataset = read.table('nome.dat',col.names=c('Ex','Ad','Tr','Gl','Op'))
write.table(dataset, file = 'dataset.txt')
  
means = tapply(kimono$value,kimono$gruppo,mean)
means = sapply(dataset,mean) #mi calcola la media per ogni col (multi dim)
                             #mean,sd,var
matr_cov = cov(dataset) #matrice covarianza
matr_cor = cor(dataset) #matrice correlazione

bartlett.test(kimono$value, kimono$gruppo)
var.test(cyto[indx_a,1],cyto[indx_b,1]) #test varianze uguali uni-dim

dnorm(0,mean = 1, sd = 2) # density function at 0 for a distribution N(1,4)
pnorm(0,1,2) # Cumulative distr: P(Z<0), with Z ~ N(1,4)
qnorm(0.95,1,2) #quantile  =z s.t. P(Z<z)=0.95, with Z ~ N(1,4)
rnorm(10,1,2) # X_1,..,X_10 ~ N(1,4) i.i.d.

# Commands to generate random numbers, obtain pdfs,
# cdf and its inverse for the most popular models:
#
# Commands rnorm(), dnorm(), pnorm(), qnorm(),
# rexp(), dexp(), pexp(), qexp(),
# runif(), dunif(), punif(), qunif(),
# rbinom(), dbinom(), pbinom(), qbinom(),
# rpois(), dpois(), ppois(), qpois(),
# rgamma(), dgamma(), pgamma(), qgamma(),
# 
# qchi(),qf(),qt() per chi-quadro,Fisher,t-stud


set.seed(123)
rep(5:6, each = 2)
rep(5:6, times = 2)
seq(2,5,len=4)
seq(2,5,by=1)
sample(x,y, replace = TRUE) #estrai x elementi da y con reimmissione = true
matrix(data = c(11,12,13,14,15,16), nrow = 2, ncol = 3, byrow = T) #riempio per righe!
solve(matrix) #inverse
solve(matrix,vect) #solve linear syst
var_gen = det(cov_matr) #compute the generalized variance
var_tot = sum(diag(cov_matr)) #compute the total variance
d2 <- mahalanobis(X, colMeans(X), cov(X)) #compute the squared mahalanobis dist
eigen(S) # $values, $vectors
exp(vect) #component-wise
sum(vect) #somma elementi del vettore 
prod(c) #prodotto elementi del vettore
round(val,2) #tronca con arrotondamento alla 2a cifra decimale
scale(dataset, center = T, scale = T) #standardizza tutte le variabili
                                      #cioè le centra(center = T) e le riscala (scale = T)

district <- c('MI',  'MI',  'VA',  'BG',  'LO', 'LO',  'CR',  'Alt', 'CR',  'MI',  
              'Alt', 'CR',  'LO',  'VA',  'MI',  'Alt', 'LO',  'MI')
district <- factor(district,levels=c('MI','LO','BG','CR','VA','Alt'))
p.adjust(pval,method = 'bonferroni') # Si ricava da pval < alfa/k --> p_corrected = k*pval
p.adjust(pval,method = 'BH')




# 2)SALVATAGGIO PLOT #############################################################

plot(rnorm(10),rnorm(10))
bmp(file = "myplot.bmp")
plot(rnorm(10),rnorm(10))
dev.off()

jpeg(file = "myplot.jpeg")
plot(rnorm(10),rnorm(10))
dev.off()

png(file = "myplot.png")
plot(rnorm(10),rnorm(10))
dev.off()

pdf(file = "myplot.pdf")
plot(rnorm(10),rnorm(10))
dev.off()

pdf(file = "myplot2.pdf", onefile = T)
plot(rnorm(10),rnorm(10))
plot(rnorm(10),rnorm(10))
plot(rnorm(10),rnorm(10))
dev.off()


# 3)GRAFICI VARI #################################################################

## 3.1)Generico ####################################################################
x11() #nuova finestra
par(mfrow = c(1,4)) #mi fa mettere 1 riga con 4 grafici sullo stessa figure

hist(campo) # deve essere campo numerico!
boxplot(campo) #deve essere numerico
boxplot(campo ~ gruppo)
boxplot(dataset[,1:2]) #boxplot delle prime due colonne
boxplot(dataset) #boxplot di tutte le colonne
boxplot(scale(x=dataset,center = T, scale=F), las=2, col='gold') #boxplot per tutte le colonne centrato!
abline() #Si aggiunge su un plot esistente!Guarda help per dettagli!
abline(h = 3) #linea orizzontale ad ascissa=3
abline(v = 5) #linea verticale ad ordinata=5
points() #come plot ma si aggiunge sopra a un plot
pairs(dataset)
matplot() #?
pie(table(variabile_factor),col=rainbow(length(levels(district)))) #fa grafico a torta di variabile cat
plot(variabile_factor) # barplot of absolute frequences
barplot(table(variabile_factor)/length(variabile_factor),ylim=c(0,1)); box() #barplot delle freq
image(x,y,z) #grafico in 2d con "heat-map" per la z
contour(x,y,z) #grafico in 2d con curve di livello per z
persp(x, y, w, theta=30, phi=30, shade=.05) #plot 3d di superficie in prospettiva settabile
lines3d(x,x, function(x,x), col='blue', lty=1) #linea 3d
biplot(pc.dataset) #plot delle proiezioni delle osservazioni sulle comp principali
text(x,y,dimnames(dataset)[[1]], cex=0.7) #mi plotta i nomi delle osservazioni al posto dei dots
rect(Bf[1,1],Bf[2,1],Bf[1,3],Bf[2,3], border='orange', lwd=2)


## 3.2)Plot ellisse3d ##############################################################
{
open3d()                    # open a new device
points3d(X, asp=1, size=4)  # plot the points
axes3d()                    # add the axes
plot3d(ellipse3d(S, centre=M, level= 0.99), alpha=0.15, add = TRUE) #alpha = intensità colore ellisse , level = proporzione di punti all'interno
}

## 3.3)Plot SimIC nelle direzioni delle componenti (ossia variabili iniziali) ######
{
x11()
k #numero di IC da plottare
#SimIC in cui ogni RIGA ha inf-center-sup
matplot(1:k,1:k,pch='',xlab='Variables',ylab='T2 for a component',
    main='Simultaneous T2 conf int for the components')
for(i in 1:k)
  segments(i,SimIC[i,1],i,SimIC[i,3],lwd=3)
points(1:k, SimIC[,2], pch=16)
# Is mu0 inside the rectangular region?
# We add it to the plot
points(1:k, mu0, lwd=3, col='orange')
}

## 3.4)Plot IC in 2dim #########################################################
{
x11()
plot(dataset, asp=1, pch=1, main='Dataset with SimIC')
ellipse(center=dataset_mean, shape=dataset_cov/n, radius=sqrt(cfr_fisher), lwd=2)
points(delta0[1], delta0[2], pch=16, col='grey35', cex=1.5) #dot for delta0
segments(SimIC[1,1],delta0[2],SimIC[1,3],delta0[2],lty=1,lwd=2,col='red')
segments(delta0[1],SimIC[2,1],delta0[1],SimIC[2,3],lty=1,lwd=2,col='red')
}

## 3.5)Worst SimIC e plot ##########################################################

#Worst direction among SimIC 
{
  worst = solve(diff_cov) %*% (diff -delta0)
  worst = worst/sqrt(sum(worst^2))
  theta_worst = atan(worst[2]/worst[1])+pi
  IC_worst #usa a = worst e calcola con codice sopra
  (IC_worst[1] < delta0%*%worst) & (delta0%*%worst < IC_worst[2]) #If true, we are inside the ellipse
  #If false, we are outside the ellipse
  x.min <- IC_worst[1]*worst
  x.max <- IC_worst[3]*worst
  m1.ort <- -worst[1]/worst[2]
  q.min.ort <- x.min[2] - m1.ort*x.min[1]
  q.max.ort <- x.max[2] - m1.ort*x.max[1]
  abline(q.min.ort, m1.ort, col='forestgreen', lty=2,lwd=1)
  abline(q.max.ort, m1.ort, col='forestgreen', lty=2,lwd=1)
  m1=worst[2]/worst[1] # worst direction
  abline(0, m1, col='grey35')
  segments(x.min[1],x.min[2],x.max[1],x.max[2],lty=1,lwd=2, col='forestgreen')
}

# 4)IPOTESI GAUSSIANITA' #########################################################

#Plot con istogramma e gaussiana teorica
{x11()
k = 1
par(mfcol = c(2,k)) #k numero di colonne da plottare una per volta
for(i in 1:k){
  hist(dataset[,i], prob=T, ylab='density', xlab=paste('Variabile', i,sep = ' '), main='Histogram of variabile')
  lines((-1000):1000 /100, dnorm((-1000):1000 /100,mean(dataset[,i]),sd(dataset[,i])), col='blue', lty=2)
  qqnorm(dataset[,i], main=paste('QQplot of V', i, sep=''))
  qqline(dataset[,i])
}
}

{qqnorm(colonna) # quantile-quantile plot (funziona su 1-dim ovviamente!)
qqline(colonna, col='red') # theoretical line
}

#qqplot per distribuzione exp (cioè guardo se campione è esponenziale)
{y <- rexp( n=1000, rate=1)
qqplot(qexp((1:1000/1000-0.5/1000)), y, col='red', xlab='Theoretical quantile', ylab='Sample Quantile', asp=1)
abline(0, 1, col='blue')
}

{#Consider the squared Mahalanobis distances of the data from the (sample) mean
#and test if they are a sample from a chi-square distribution
  
#Recall:
    # Theorem: if X~N(mu,Sigma) r.v. in R^p, det(Sigma)>0
    #          then d2(X,mu)=(X-mu)'Sigma^-1(X-mu) ~ Chi-sq(p)
d2 = mahalanobis(dataset, colMeans(dataset), cov(dataset))
n = dim(dataset)[1] #numero di osservazioni 
x11(width=13)
par(mfrow=c(1,2))
hist(d2, prob=T, main='Histogram of the Mahalanobis dist.',xlab='d2',ylab='density', col='grey84')
lines(0:2000/100, dchisq(0:2000/100,numero_colonne), col='blue', lty=2, lwd=2)
qqplot(qchisq((1:n - 0.5)/n, df = numero_colonne), d2, main='QQplot of (sample) d2',xlab='theoretical quantiles Chi-sq(2)',
     ylab='sample quantiles')
abline(0, 1)
 
d2.class <- cut(d2, qchisq((0:10)/10, df = numero_colonne))
d2.freq <- table(d2.class)
chisq.test(x = d2.freq, p = rep(1/10, 10), simulate.p.value = T)  #H0: è chi-sq 
}

shapiro.test(dataset$campo)$p

mcshapiro.test(dataset_multi_dim)

{#mcshapiro senza outliers rispetto la d2 dist
  d2 = mahalanobis(dataset, colMeans(dataset), cov(dataset))
  dataset_senza_outliers = dataset[which(d2<7.5),]
  mcshapiro.test(dataset_senza_outliers)
  
}

{#Box-Cox 1-dim
  #box_cox <- function(x,lambda)
  #{
  #  if(lambda!=0)
  #    return((x^lambda-1)/lambda)
  #  return(log(x))
  #  1
  #}
  # lambda<1: observations <1 are "spread", observations >1 are "shrinked"
  lambda_opt = powerTransform(colonna)
  colonna_bc = bcPower(colonna, lambda_opt$lambda) 
}

{#Box-Cox multi-dim
  # lambda<1: observations <1 are "spread", observations >1 are "shrinked"
  lambda_opt_vect = powerTransform(dataset) #cbind(col1,col2) per due colonne
  dataset_bc = bcPower(dataset, lambda_opt_vect$lambda) 
}



bartlett.test(kimono$value, kimono$gruppo) #varianze uguali nei gruppi


#5)PCA (lab 3) ##################################################################

#NB: Start with visualization : boxplot, pairs, image(cov_matr) 

#Plot con ellisse e direzioni principali per 2dim
x11()
plot(dataset_2dim, asp=1, xlab='Var 1', ylab='Var 2',pch=20)
ellipse(mean, S, radius, add=T,lwd=3, col='red') #S sarebbe la matr cov
                                                 #al posto di radius metti level = 0.55 per mettere la proporzione(0.55) di punti all'interno
abline(a = M[2] - eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]*M[1], b = eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1], lty = 2, col = 'dark red', lwd = 2)
abline(a = M[2] - eigen(S)$vectors[2,2]/eigen(S)$vectors[1,2]*M[1], b = eigen(S)$vectors[2,2]/eigen(S)$vectors[1,2], lty = 2, col = 'red', lwd = 2)

#Comandi pca
pc.dataset = princomp(dataset, scores=T) #$sd stand dev of components
summary(pc.dataset)

proportion_variance_per_comp = pc.dataset$sdev^2/sum(pc.dataset$sdev^2)
comulat_proportion_varaince = cumsum(pc.dataset$sdev^2)/sum(pc.dataset$sdev^2)
loads = pc.dataset$loading #empty are zeros or near
scores = pc.dataset$scores  #proiezioni osservazioni su componenti principali
pc_i = loads[i,] #vettore della comp princ i-esima
new_scores = loads%*%new_obs

#Grafico dei loads delle prime 8 comp princip
{x11() 
par(mfcol = c(4,2))
for(i in 1:8)
  #taglio i loads sotto una soglia minima
  barplot(ifelse(abs(loads[,i]) < soglia_minima, 0, loads[,i]), ylim = c(-1, 1), main=paste("Loadings PC",i, sep = ' '))
}
#Plot per explained variance
{x11()
layout(matrix(c(2,3,1,3),2,byrow=T))
plot(pc.dataset, las=2, main='Principal components')
barplot(sapply(dataset,sd)^2, las=2, main='Original Variables', ylab='Variances') #riga 194 lab3
plot(cumsum(pc.dataset$sd^2)/sum(pc.dataset$sd^2), type='b', axes=F, xlab='number of components',
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(dataset),labels=1:ncol(dataset),las=2)
}

#Plot scores : plot delle proiezioni sulle prime 2 comp princ
plot(scores[,1:2]) 
biplot(pc.dataset) #analogo ma più fiko??

#Confronto dataset e proiezioni su comp princ
boxplot(dataset, las=2, col='gold', main='Original variables')
boxplot(data.frame(scores.dataset), las=2, col='gold', main='Principal components')

#Standardizzazione variabili (ripeti analisi fatte sopra!)
dataset.std <- scale(dataset)
dataset.std <- data.frame(dataset.std)



#Riga 415 lab3 per matplot di "Projection on the space generated by the k-th principal component"
#Riga 429 lab3 per matplot di "Projection on the space generated by the first k principal components"



#6)TEST/IC MEDIA E VARIANZA DI POP GAUSSIANA ###################################


tstat <- (sample.mean - mean.H0)/(sample.sd/sqrt(n))
cfr.t <- qt( 1 - alpha/2, n-1 ) #NB: cambiare i df (ossia n-1) della S se necessario!
abs(tstat) < cfr.t
pval <- ifelse(tstat >= 0, (1 - pt(tstat, n-1))*2, pt(tstat, n-1)*2)
IC <- c(inf = sample.mean - sample.sd/sqrt(n) * qt( 1 - alpha/2, n-1 ),
        center = sample.mean,
        sup = sample.mean + sample.sd/sqrt(n) * qt( 1 - alpha/2, n-1 ))



# 7)TEST MEDIA #################################################################
##7.1)Teoria  ##################################################################

#Hotelling:
# X~Np(mu,??) , det(??)>0
# W =somma(1:m)of(Zi'Zi)~Wish(??/m,m)
# X indip W
# Then:
# (m-p+1)/(m*p)*(sample_mean_X-mu)'W^-1(sample_mean_X-mu)~F(p,m-p+1)

#In generale:
# p = covariate di X cioè dim(dataset)[2]
# n = osservazioni di X cioè dim(dataset)[1]
# m = degree of freedom di W^-1
# W = quasi sempre è stimatore della matr_cov di sample_mean_X
#     QUINDI SARA' (S/n)! Occhio alle Spool in cui n = (1/n1+1/n2)

#Regione di confindenza:
# - centro in dataset_mean
# - direzioni assi principali sono gli autovett della matrice W^-1 cioè stimatore della
#   matr_cov di ciò che sto considerando (ex: per X sarà S, per X_sample_mean sarà S/n)
#   quindi eigen(dataset_cov/n)$vectors (nel caso di X_sample_mean)
# - raggio è sqrt(cfr_fisher) o sqrt(cfr_chisq) dai casi
# - lunghezza semiassi = raggio*sqrt(eigen(dataset_cov/n)$values) (nel caso di "stimatore di S" = cov/n)

#-If the sample mean (red dot) falls in the rejection region ( region OUTSIDE blue ellipse) 
#  we reject H0, if it falls in the acceptance region (inside blue ellipse) we accept H0
#-If the m0 (blue dot) falls outside the conf region (outside red ellipse) we reject H0,
#  if it falls in the confidence region (inside red ellipse) we accept H0


#Stimatori di ??:
# -X~Np(mu,??) ==> S , df = n-1
# -X1,..Xn~Np(mu,??), X_sample_mean ==> S/n , df = n-1 
# -differenza tra medie di gruppi con n1!=n2 ==> Spool/(1/n1+1/n2) , df =(n1-1)+(n2-1)
# -da finire


## 7.2)Test e ic media ####

###7.2.1)Pop gaussiana ##############################################################
{
  #Test media = mu0
  alpha #da inserire
  k #da inserire numero di test di Bonferroni
  mu0 #ipotesi, da inserire
  n = dim(dataset)[1] #numero osservazioni, da inserire
  p = dim(dataset)[2] #covariate
  df = n-1 # = m cioè degree of freedom di S (matr_cov),nel caso "classico" sono n-1
  dataset_mean = sapply(dataset,mean)
  dataset_cov = sapply(dataset,cov)
  inv_cov = solve(dataset_cov)
  T2 = n*(dataset_mean-mu0) %*% inv_cov %*% (dataset_mean-mu0)
  cfr_fisher = ((df)*p/(df-p+1))*qf(1-alpha,p,df-p+1) #NB: se si sostituisce df=n-1, si ottiene caso classico (coroll sotto Hotelling nelle lezioni)
                                                    # Oss: cfr_fisher = radius^2 dell'ellisse
  
  T2 < cfr_fisher #if true, accetto H0!
  pval = 1-pf(T2*(df-p+1)/((df)*p), p, df-p+1)
  
 #SimIc su direzione particolare
  a #da definire vettore su cui proiettare
  SimIC = cbind(inf = a%*%dataset_mean - sqrt(cfr_fisher*(a%*%dataset_cov%*%a)/n),
                center = a%*%dataset_mean,
                sup = a%*%dataset_mean + sqrt(cfr_fisher*(a%*%dataset_cov%*%a)/n))
  #SimIc su direzioni degli assi 
  SimIC_assi = cbind(inf = dataset_mean - sqrt(cfr_fisher*diag(dataset_cov)/n),
                     center = dataset_mean,
                     sup = dataset_mean + sqrt(cfr_fisher*diag(dataset_cov)/n))
  
  
  a #da definire vettore su cui proiettare
  k #da inserire numero di test di Bonferroni
  cfr_t = qt(1-alpha/(2*k),df)
  BonfIC = cbind(inf = a%*%dataset_mean - cfr_t*sqrt( (a%*%dataset_cov%*%a)/n),
                 center = dataset_mean,
                 sup = a%*%dataset_mean + cfr_t*sqrt((a%*%dataset_cov%*%a)/n))
}  
  
#### 7.2.1.1)Plot cfr SOLO per p=2 !!!!!!!! #################################################
  x11()
  #Rejection region = outside blue ellipse, centred in mu0 (blue dot)
  plot(dataset, asp = 1)
  ellipse(mu0, shape=dataset_cov/n, sqrt(cfr_fisher), col = 'blue', lty = 2, center.pch = 16)
  points(dataset_mean[1], dataset_mean[2], pch = 16, col ='red', cex = 1.5)# Red point for sample mean
  #Conf region = inside red ellipse, centred in the sample mean (red dot)
  ellipse(dataset_mean, dataset_cov/n, sqrt(cfr_fisher), col = 'red', lty = 2, lwd=2, center.cex=1)
  
  #-If the sample mean (red dot) falls in the rejection region ( region OUTSIDE blue ellipse) 
  #  we reject H0, if it falls in the acceptance region (inside blue ellipse) we accept H0
  #-If the m0 (blue dot) falls outside the conf region (outside red ellipse) we reject H0,
  #  if it falls in the confidence region (inside red ellipse) we accept H0
  
 

###7.2.2) Asintotici #################################################################
{
  mu0 #ipotesi, da inserire
  n = dim(dataset)[1] #numero osservazioni, da inserire
  p = dim(dataset)[2] #covariate
  dataset_mean = sapply(dataset,mean)
  dataset_cov = sapply(dataset,cov)
  inv_cov = solve(dataset_cov)
  T2 = n*(dataset_mean-mu0) %*% inv_cov %*% (dataset_mean-mu0)
  cfr_chisq = qchisq(1-alpha, p)
  
  T2 < cfr_chisq
  pval = 1-pchisq(T2,p)
  
#### 7.2.2.1)Plot cfr SOLO per p=2 !!!!!!!! #################################################
  x11()
  #Rejection region = outside blue ellipse, centred in mu0 (blue dot)
  plot(dataset, asp = 1)
  ellipse(mu0, shape=dataset_cov/n, sqrt(cfr_chisq), col = 'blue', lty = 2, center.pch = 16)
  points(dataset_mean[1], dataset_mean[2], pch = 16, col ='red', cex = 1.5)# Red point for sample mean
  #Conf region = inside red ellipse, centred in the sample mean (red dot)
  ellipse(dataset_mean, dataset_cov/n, sqrt(cfr_chisq), col = 'red', lty = 2, lwd=2, center.cex=1)
}


### 7.2.3) Mono-dir con gauss  ########################################################
{
  #Test in una singola direzione per pop multidim gaussiana
  n = dim(dataset)[1] #numero di osservazioni
  dataset_mean = sapply(dataset,mean)
  dataset_cov = sapply(dataset,cov)
  df = n-1 #degree of freedom di stimatore di varianza della sample mean
  delta0 #da definire
  alpha #da definire
  a #da definire, UNICA direzione che considero
  T0 = ((a%*%dataset_mean)-delta0)/sqrt(a%*%dataset_cov%*%a)*sqrt(n)

  
  #H0 : a*mu <= delta0 --> cfr_t = qt(1-alpha,df)
  #                        We reject if T0 > cfr_t
  #                        pval = 1- pt(T0, df)
  
  #HO : a*mu > delta0 -->  cfr_t = qt(1-alpha,df)
  #                        We reject if T0 < cfr_t
  #                        pval = pt(T0, df)
  
  #H0 : a*mu = delta0 -->  cfr_t = qt(1-alpha/2,df)
  #                        We reject if abs(T0) > cfr_t
  #                        pval = ifelse(T0 >= 0, (1 - pt(T0,df))*2, pt(T0,df)*2)
  
  
}



# 8) PAIRED DATA lab6##################################################################

#n = numero di unità
#p = numero di features
#q = numero di misure ripetute = 2

## 8.1) Ipotesi #####################################################################
# Calcola le DIFFERENZE e devono essere Np(delta,sigmadD) 

## 8.2) Test e ic differenza medie ####################################################### 
  
#NB: ATTENZIONE A DEFINIRE n,p,q !!!!!

delta0=c(0,0)
alpha #da definire
n 
p
q
df = n-1 #df di stimatore della varianza della media (classico è S/n,df=n-1)
diff_mean = sapply(diff,mean)
diff_cov = sapply(diff,cov)
inv_cov = solve(diff_cov)

T0 = n * (diff_mean - delta0)%*%inv_cov%*%(diff_mean- delta0)
cfr_fisher = ((df*p)/(df-p+1))*qf(1-alpha,p,df-p+1)

T0 < cfr_fisher #Reject if false
pval = 1 - pf(T0*(df-p+1)/(df*p),p,df-p+1)

#SimIC su direzione particolare
{
a #da definire vettore su cui proiettare
SimIC =cbind(inf =  a%*%diff_mean -sqrt(cfr.fisher* (a%*%diff_cov%*%a)/n),
             center = a%*%diff_mean,
             sup =  a%*%diff_mean +sqrt(cfr.fisher* (a%*%diff_cov%*%a)/n) )
}

#SimIC nelle direzioni degli assi
SimIC =cbind(inf =  diff_mean -sqrt(cfr.fisher* diag(diff_cov)/n),
             center = diff_mean,
             sup =  diff_mean -sqrt(cfr.fisher* diag(diff_cov)/n) )

#Worst direction among SimIC 
{
worst = solve(diff_cov) %*% (diff -delta0)
worst = worst/sqrt(sum(worst^2))
theta_worst = atan(worst[2]/worst[1])+pi
IC_worst #usa a = worst e calcola con codice sopra
(IC_worst[1] < delta0%*%worst) & (delta0%*%worst < IC_worst[2]) #If true, we are inside the ellipse
                                                                #If false, we are outside the ellipse
x.min <- IC_worst[1]*worst
x.max <- IC_worst[3]*worst
m1.ort <- -worst[1]/worst[2]
q.min.ort <- x.min[2] - m1.ort*x.min[1]
q.max.ort <- x.max[2] - m1.ort*x.max[1]
abline(q.min.ort, m1.ort, col='forestgreen', lty=2,lwd=1)
abline(q.max.ort, m1.ort, col='forestgreen', lty=2,lwd=1)
m1=worst[2]/worst[1] # worst direction
abline(0, m1, col='grey35')
segments(x.min[1],x.min[2],x.max[1],x.max[2],lty=1,lwd=2, col='forestgreen')
}

#Bonferroni ic in direzione generica
{
  a #da definire vettore su cui proiettare
  k #da inserire numero di test di Bonferroni
  df = n-1 
  cfr_t = qt(1-alpha/(2*k),df)
  BonfIC = cbind(inf = a%*%diff_mean - cfr_t*sqrt((a%*%dataset_cov%*%a)/n),
                 center = a%*%diff_mean,
                 sup = a%*%diff_mean + cfr_t*sqrt((a%*%dataset_cov%*%a)/n))
}

#Bonferroni ic sulle componenti
{
  k = p #da inserire numero di test di Bonferroni (default guardo TUTTE le componenti)
  df = n-1 
  cfr_t = qt(1-alpha/(2*k),df)
  BonfIC = cbind(inf = diff_mean - cfr_t*sqrt(diag(dataset_cov)/n),
                 center = diff_mean,
                 sup = diff_mean + cfr_t*sqrt(diag(dataset_cov)/n))
}

## 8.3) Plot cfr 2dim ###############################################################
{
#Conf region : inside the ellipse
x11()
plot(diff, asp=1, pch=1, main='Dataset of the Differences')
ellipse(center=diff_mean, shape=diff_cov/n, radius=sqrt(cfr_fisher), lwd=2)
points(delta0[1], delta0[2], pch=16, col='grey35', cex=1.5) #dot for delta0
segments(SimIC[1,1],delta0[2],SimIC[1,3],delta0[2],lty=1,lwd=2,col='red')
segments(delta0[1],SimIC[2,1],delta0[1],SimIC[2,3],lty=1,lwd=2,col='red')
}

# 9) REPEATED MEASURES with p=1 q>=2 lab6 #########################################

## 9.1)Ipotesi #####################################################################

#NB: ATTENZIONE A DEFINIRE n (num osservaz),q(ripetizioni delle misurazioni),p(numero covariate)
#df = n-1 (=m in Hotelling), q-1 (=p in Hotelling) 

#Ipotesi : (X1i,..Xqi) ~ Nq(mu, matr_sigma) quindi mcshapiro su dataset nxq



## 9.2) Test e ic ###################################################################

alpha # da definire
delta0 #a definire
C #contrast matrix (q-1)xq, righe sono le differenze (usa matrix(...,byrow=T))
m_c = C %*% dataset_mean
cov_c = C %*% dataset_cov %*% t(C)
invcov_c = solve(cov_c)

T0 = n* t( m_c - delta0 ) %*% invcov_c %*% ( m_c - delta0 )
cfr_fisher = (n-1)*(q-1)/(n-q+1)*qf(1-alpha,q-1,n-q+1)

T0 < cfr_fisher #If false,we reject
pval = 1 - pf(T0*(n-q+1)/( (q-1)*(n-1) ),q-1,n-q+1)

#SimIC nella direzione degli assi
{
  SimIC =cbind(inf =  m_c -sqrt(cfr.fisher* diag(cov_c)/n),
               center = m_c,
               sup =  m_c +sqrt(cfr.fisher* diag(cov_c)/n) )
}


#Bonferroni IC nella direzione degli assi
{
 k = q-1 #numero di righe della contrast matrix
 cfr_t = qt(1-alpha/(2*k),n-1)
 BonfIC = cbind(inf = m_c - cfr_t*sqrt(diag(cov_c)/n),
                center = m_c,
                sup = m_c + cfr_t*sqrt(diag(cov_c)/n))
}

# 10)TWO INDEP GAUSS POP lab6#######################################################

## 10.1) Ipotesi #####################################################################
# Spool df = n1+n2-2
# Ipotesi : verifica che le due popolazioni siano gaussiane e con uguale varianza!
#           L'indipendenza la da il problema(?)
bartlett.test(kimono$value, kimono$gruppo)
mcshapiro.test(dataset_multi_dim)



## 10.2) Test e ic ###################################################################
n1 = dim(dataset1)[1]
n2 = dim(dataset2)[1]
p = dim(dataset1)[2] #(p1 = p2 !)

dataset1_mean = sapply(dataset1,mean)
dataset2_mean = sapply(dataset2,mean)
dataset1_cov = sapply(dataset1,cov) #se non funziona, usa cov(dataset1)
dataset2_cov = sapply(dataset2,cov)
Sp = ( (n1-1)*dataset1_cov + (n2-1)*dataset2_cov)/(n1+n2-2)

alpha #da definire
delta0 #da definire
Spinv = solve(Sp)
T0 = n1*n2/(n1+n2) * (dataset1_mean - dataset2_mean - delta0) %*% Spinv %*% (dataset1_mean - dataset2_mean - delta0)
cfr_fisher = (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)

T0 < cfr_fisher #If false, we reject
pval = 1 - pf(T0/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)

#SimIC nella direzione degli assi
{
  SimIC =cbind(inf =  dataset1_mean - dataset2_mean -sqrt(cfr_fisher* diag(Sp)*(1/n1+1/n2)),
               center = dataset1_mean - dataset2_mean,
               sup =  dataset1_mean - dataset2_mean +sqrt(cfr_fisher* diag(Sp)*(1/n1+1/n2)) )
}

#BonfIC nella direzione degli assi
{
  k = p #numero di test = p se guardo gli assi
  cfr_t = qt(1-alpha/(2*k),n1+n2-2)
  BonfIC = cbind(inf = dataset1_mean - dataset2_mean - cfr_t*sqrt(diag(Sp)*(1/n1+1/n2)),
                 center = dataset1_mean - dataset2_mean,
                 sup = dataset1_mean - dataset2_mean + cfr_t*sqrt(diag(Sp)*(1/n1+1/n2)) )
}


# 11) PROPORTIONS/BERNULLI penultimo es lab6 #########################################################
#era un es sulle indep gauss populations ma meglio tenerlo come esempio specifico

## 11.1) Ipotesi #####################################################################
# Ipotesi : non ha controllato un cazzo di niente...usa tcl (?)
## 11.2)Test e ic ###################################################################

n1 <- dim(allergy)[1]
n2 <- dim(noallergy)[1]
p <- dim(noallergy)[2]
x.mean1 <- sapply(allergy, mean)
x.mean2 <- sapply(noallergy, mean)
p.hat <- (x.mean1*n1+x.mean2*n2)/(n1+n2)
x.var <- (p.hat*(1-p.hat))
# Test: H0.i: mu.i1 == mu.i2 vs H1.i: mu.i1 != mu.i2
z.i <- (x.mean1-x.mean2)/sqrt(x.var*(1/n1+1/n2))
p.i <- ifelse(z.i<0, 2*pnorm(z.i),2*(1-pnorm(z.i)))
#Bonferroni test
k <- 520
which(p.i*k<.01)
# or
p.Bf <- p.adjust(p.i, method='bonferroni')
which(p.Bf<.01)
# Benjamini-Hockberg (control the false discovery rate)
p.BH <- p.adjust(p.i, method='BH')
which(p.BH<.01)

# 12)ANOVA ONE-WAY lab7#################################################################
# p = 1 , g>=1

## 12.1) Ipotesi #####################################################################
n = dim(dataset)[1] # total number of obs.
ng <- table(dataset$group) # number of obs. in each group
gruppi <- levels(dataset$group) # levels of the treatment
g <- length(gruppi) # number of levels (i.e., of groups)

plot(dataset$value,dataset$gruppo, xlab='treat', col='grey85')
shapiro.test(dataset$value[ dataset$gruppo==gruppi[1] ])$p #per gruppo 1 , ripeti per tutti i gruppi
var(dataset$value[ dataset$gruppo==gruppi[1] ]) #per gruppo 1, ripeti per tutti i gruppi
bartlett.test(dataset$value, dataset$gruppo)

SSres = sum(residuals(fit)^2)
df_res = fit$df.residual #n-g
S = SSres/df_res
media_g = tapply(dataset$value,dataset$gruppo,mean)

## 12.2) Analisi #####################################################################

fit = aov(dataset$value ~ dataset$gruppo )
summary(fit)  
#How to read the summary:
  # Df Sum Sq Mean Sq F value Pr(>F)
  # treat (g-1) SStreat SStreat/(g-1) Fstatistic p-value [H0: tau.i=0 for every i]
  # Residuals (n-g) SSres SSres/(n-g)
  
#If we reject the test, i.e., we have evidence to state that the
# treatment (feed supplement) has an effect on the growth rate
# of chicken , we want to understand which supplement is responsible for this.
#To see this, we need to do g*(g-1)/2 comparisons.


## 12.3)Confronti ###################################################################
#SStot = SSmean + SStreat + SSres
#Df: n =    1   +  g-1    + n-g 

#SScentered = SStreat + SSres
#Df:  n-1   =   g-1   + n-g


{# BonfIC Differenza medie di due gruppi :

media_g = tapply(dataset$value,dataset$gruppo,mean)
ng <- table(dataset$group)
SSres = sum(residuals(fit)^2)
df_res = fit$df.residual #n-g
S = SSres/df_res

ic_bonf = cbind(inf = media_g[1]- media_g[2] - qt(1-alpha/(2*k), n-g) * sqrt(S * (1/ng[1] + 1/ng[2]))  ,
                center = media_g[1]- media_g[2],
                sup = media_g[1]- media_g[2] + qt(1-alpha/(2*k), n-g) * sqrt(S * (1/ng[1] + 1/ng[2])))
}

{#BonfIC di tutte le differenze tra i gruppi (g*(g-1)/2 confronti)
Mediag = tapply(dataset$value,dataset$gruppo,mean)
k = g*(g-1)/2
alpha #da definire
ICrange=NULL
gruppi = levels(dataset$group)
ng = table(dataset$group) 
for(i in 1:(g-1)) {
  for(j in (i+1):g) {
      print(paste(gruppi[i],"-",gruppi[j]))
      print(as.numeric(c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )) ,
                         Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] ))) ))
      ICrange=rbind(ICrange,as.numeric(c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                                         Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
  }
}
x11(width = 14, height = 7)
par(mfrow=c(1,2))
plot(feed, weight, xlab='treat', ylab='weight', col = rainbow(6), las=2)
h <- 1
plot(c(1,g*(g-1)/2),range(ICrange), pch='',xlab='pairs treat')
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
}


{
#Ic_bonf media di un gruppo:
k #da definire
S_i = var(dataset$value[dataset$gruppo == gruppi[i]])
ic_bonf = cbind(inf = media_g[i] - qt(1-alpha/(2*k),ng[i]-1) * sqrt(S_i/ng[i]), 
                center = media_g[i],
                sup = media_g[i] + qt(1-alpha/(2*k),ng[i]-1) * sqrt(S_i/ng[i]))
}

{
#Ic_bonf varianza di un gruppo:
ic_bonf = cbind(inf = (ng[i]-1)*S_i/qchisq(1-alpha/(2*k),ng[i]-1),
                center = 0 ,
                sup = (ng[i]-1)*S_i/qchisq(alpha/(2*k),ng[i]-1) )
}



# 13)MANOVA ONE-WAY lab7 ###############################################################
#n osservazioni, g>=2 gruppi, p>=2 features

## 13.1) Ipotesi #####################################################################

#Ipotesi 1: normalità (multivariata) in ogni gruppo
gruppi = factor(dataset$gruppo, labels=c(unique(dataset$gruppo)))
i1 = which(dataset$gruppo == gruppi[1]) #estraggo indici gruppo 1
mcshapiro.test(dataset[i1])
n1 <- length(i1) #fare anche per gli altri gruppi
n <- n1+n2+n3

#Ipotesi 2: same covariance structure
S = cov(dataset)
S1 <- cov(dataset[i1,])
S2 <- cov(dataset[i2,])
S3 <- cov(dataset[i3,])

bartlett.test(dataset$value, dataset$gruppo) #usalo sul modello ridotto!!

# Note: We can verify the assumptions a posteriori on the residuals of
# the estimated model



## 13.2) Analisi #####################################################################

fit = manova(as.matrix(dataset) ~ gruppi) #gli passo gruppi e NON dataset$gruppo
summary.manova(fit,test="Wilks")

# If Reject the test, i.e., we have statistical evidence to state that
# the factor "gruppo" has an effect on the mean features
# of the dataset, we want to understand who is the responsible for this.

#Via ANOVA: for each of the p variables we perform an ANOVA test
# to verify if the membership to a group has influence
# on the mean of the variable (we explore separately the
# p axes directions in R^p)
summary.aov(fit) # guardiamo p modelli anova in cui ho "dataset$variabile ~ dataset$gruppo"

# Note: this analysis does NOT say:
# a) which group differ
# b) with respect to which variables the groups in (a) differ
# => As for the ANOVA, we build confidence intervals (many more!)



## 13.3)Confronti ###################################################################

{#BonfIC (nel codice si è supposto p = 4 !)
  alpha <- 0.05
  k <- p*g*(g-1)/2
  qT <- qt(1-alpha/(2*k), n-g)
  W <- summary.manova(fit)$SS$Residuals
  m <- sapply(dataset,mean) # estimates mu
  m1 <- sapply(dataset[i1,],mean) # estimates mu.1=mu+tau.1
  m2 <- sapply(dataset[i2,],mean) # estimates mu.2=mu+tau.2
  m3 <- sapply(dataset[i3,],mean) # estimates mu.3=mu+tau.3
  
  inf12 <- m1-m2 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
  sup12 <- m1-m2 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
  inf13 <- m1-m3 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n3) )
  sup13 <- m1-m3 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n3) )
  inf23 <- m2-m3 - qT * sqrt( diag(W)/(n-g) * (1/n2+1/n3) )
  sup23 <- m2-m3 + qT * sqrt( diag(W)/(n-g) * (1/n2+1/n3) )
  
  CI <- list(setosa_versicolor=cbind(inf12, sup12), setosa_virginica=cbind(inf13, sup13), versicolor_virginica=cbind(inf23, sup23))
  CI
  
  #lab7 riga 377 per plot molto lunghi o pdf pag 17/18
}


# 14) TWO-WAY ANOVA lab7  ##########################################################
# p=1 , g>=2 , b>=2

## 14.1) Ipotesi #####################################################################

#Controllare alla fine sul modello finale!
#Di solito rimane un modello ridotto del tipo anova e si controlla gaussianità e 
# matr cov
shapiro.test(km[ benz==levels(benz)[1] ])$p
bartlett.test(dataset$value, dataset$gruppo)

## 14.2) Pre-processing ##############################################################

dataset$gruppo = factor(dataset$gruppo) #per treat1 e treat2

## 14.3) Analisi #####################################################################


{#We start from the complete model
fit.aov2.int <- aov(km ~ distr + benz + distr:benz,dataset)
summary.aov(fit.aov2.int)

# Test:
# 1) H0: gamma.11 = gamma.12 = gamma.21 = gamma.22 = 0 vs H1: (H0)^c i.e.
#    H0: There is no significant interaction between the factors station
#        and gasoline in terms of performances
#    H1: There exists a significant interaction between the factors station
#        and gasoline in terms of performances
#
# 2) H0: tau.1 = tau.2 = 0 vs H1: (H0)^c i.e.
#    H0: The effect "gas station" doesn't significantly influence performances
#    H1: The effect "gas station" significantly influences performances
#
# 3) H0: beta.1 = beta.2 = 0 vs H1: (H0)^c i.e.,
#    H0: The effect "gasoline" doesn't significantly influence performances
#    H1: The effect "gasoline" significantly influences performances

#Si toglie l'interazione se non è significativa => i df_interact finiscono nei residui!

# Additive model:
# X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N(0,sigma^2),
# i=1,2 (effect station), j=1,2 (effect gasoline)
fit.aov2.ad <- aov(km ~ distr + benz, dataset)
summary.aov(fit.aov2.ad)

#Se c'è qualcosa di non significativo, si toglie => i df_treat che si toglie finiscononei residui!

fit.aov.reduced = aov(km ~ distr ,dataset )
summary(fit.aov.reduced)
}

{#Global test per significatività dei due treatments (modello senza interazione)
  #Tratto i due treatments come un unico treat : sommo i SS e i df!
  SSdistr <- sum(n*b*(Mdistr - M)^2) # or from the summary from the "Mean Sq" column
  SSbenz <- sum(n*g*(Mbenz - M)^2) # or from the summary from the "Mean Sq" column
  SSres <- sum(fit.aov2.ad$residuals^2) # or from the summary from the "Mean Sq" column
  Ftot <- ( (SSdistr + SSbenz) / ((g-1)+(b-1)))/(SSres / (n*g*b-g-b+1))
  Ptot <- 1 - pf(Ftot, (g-1)+(b-1), n*g*b-g-b+1) # attention to the dof!
  Ptot
}


## 14.4)Confronti ###################################################################

{#IC per diff medie sul reduced add model in cui rimane solo treat b
  #1)Occhio ad aggiustare il fattore finale se i gruppi hanno numerosità diverse!
  #2)Se sono bonf int, aggiusta k
  k = 1
  alpha #da definire
  SSres = sum(fit.aov.reduced$res^2)
  df = (n*g-1)*b # se è rimasto treat g, aggiustare !!!!!!!
  IC <- cbind(inf = Mbenz[1] - Mbenz[2] - qt(1-alpha/(2*k), df) * sqrt(SSres/(df) *(1/(n*g) + 1/(n*g))),
              center = Mbenz[1] - Mbenz[2],
              mean = Mbenz[1] - Mbenz[2] + qt(1-alpha/(2*k), df) * sqrt(SSres/(df) *(1/(n*g) + 1/(n*g))))
}

{#BonfIC su stime media e varianza nei gruppi del reduced model
  df <- fit3$df #df dei residui
  Spooled <- sum(fit3$res^2)/df
  means <- as.vector(tapply(penguins$weight, penguins$species, mean))
  names(means) <- levels(species)
  n #numerosità nei gruppi del treat che considero. Aggiustare nell'IC se sono diverse tra i gruppi!!!
  alpha <- 0.10
  k <- 4 # g + 1 = 4 (g Conf Int for the means and 1 for the variance which is equal in the groups for hyp!)
  BF <- rbind(cbind(means - sqrt(Spooled / n) * qt(1 - alpha / (2*k), df),
                    means + sqrt(Spooled / n) * qt(1 - alpha / (2*k), df)),
              c(Spooled * df / qchisq(1 - alpha / (2*k), df),
                Spooled * df / qchisq(alpha / (2*k), df)))
  rownames(BF)[4] <- 'Var.'
  BF
  
  
  
}





# 15) TWO-WAY MANOVA ###############################################################
# p>=1, b>=2 , g>=2

##15.1) Pre-processing ##############################################################
Ex <- factor(plastic$Ex, labels=c('L','H')) # Treat.1
Ad <- factor(plastic$Ad, labels=c('L','H')) # Treat.2

ExAd <- Ex
levels(ExAd) <- c('LL','LH','HL','HH')
ExAd[Ex=='L' & Ad=='L'] <- 'LL'
ExAd[Ex=='L' & Ad=='H'] <- 'LH'
ExAd[Ex=='H' & Ad=='L'] <- 'HL'
ExAd[Ex=='H' & Ad=='H'] <- 'HH'

{#Grafico esplorativo su prima variabile , da ripetere per tutte le p variabili
boxplot(plastic3[,1]~ExAd, main='Model with Interac. Extrusion+Additive (Tear Resistance)', ylab='Tr', col='grey95')
boxplot(plastic3[,1]~Ex,   main='Only Factor Extrusion'  , ylab='Tr', col=c('red','blue'))
}
##15.2) Ipotei ######################################################################

#Controllare gaussianità multi-dim e omogeneity matr covarianze (qualitativamente)

mcshapiro.test(plastic3[ ExAd==levels(ExAd)[1], ])$p #Ripetere per tutte le combinaz di gruppi
mcshapiro.test(plastic3[ ExAd==levels(ExAd)[2], ])$p 

S1 <- cov(plastic3[ ExAd==levels(ExAd)[1], ])
S2 <- cov(plastic3[ ExAd==levels(ExAd)[2], ])
S3 <- cov(plastic3[ ExAd==levels(ExAd)[3], ])
S4 <- cov(plastic3[ ExAd==levels(ExAd)[4], ])
x11(width=21)
par(mfrow=c(1,4))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3,S4), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3,S4), (0:100)/100, na.rm=TRUE))
image(S3, col=heat.colors(100),main='Cov. S3', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3,S4), (0:100)/100, na.rm=TRUE))
image(S4, col=heat.colors(100),main='Cov. S4', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3,S4), (0:100)/100, na.rm=TRUE))

## 15.3) Analisi #####################################################################

### Model with interaction (complete model):
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N_p(0,Sigma), [p=3]
### i=1,2 (effect Extrusion), j=1,2 (effect Additive),
### X.ijs, mu, tau.i, beta.j, gamma.ij in R^3
fit <- manova( as.matrix(plastic3) ~ Ex + Ad + Ex:Ad)
summary.manova(fit, test="Wilks")


### Model without interaction (additive model):
### X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N_p(0,Sigma), [p=3]
### i=1,2 (effect Extrusion), j=1,2 (effect additive),
### X.ijs, mu, tau.i, beta.j, in R^3
fit2<- manova( as.matrix(plastic3) ~ Ex + Ad)
summary.manova(fit2, test="Wilks")

#Per guardare QUALI treats hanno effetto sulle variabili, guardiamo il summary
# dell'anova associato all'ultimo modello manova trovato cioè:

summary.aov(fit2) #per ogni var vedo quale treat è influente
#Trovati quali treats hanno effetto tra i due (treat1 e treat2), guardo con dei
# BonfIC quale effetto hanno i g tipi di treat1 e i b tipi di treat2 sulle n variabili

## 15.4) Confronti ###################################################################

{#BonfIC per effetto
  alpha #da definire
  g 
  b
  p
  n #=numerosità in ciascun gruppo : aggiustare se i gruppi non sono equlibrati!!!!
  N = n*g*b # = totale osservazioni , aggiustare se i gruppi non sono equilibrati!!!
  W = summary.manova(fit2)$SS$Residuals
  k = g*(g-1)/2*p + b*(b-1)/2*p
  df = g*b*n-g-b+1 #aggiustare in base al modello che si usa! Default : additive
  qT <- qt(1 - alpha / (2 * k), df)
  
  mExL  <- sapply(plastic3[Ex=='L',],mean)
  mExH  <- sapply(plastic3[Ex=='H',],mean)
  infEx <- mExH-mExL - qT * sqrt( diag(W)/(df) * (1/n+1/n) )
  centerEX = mExH-mExL
  supEx <- mExH-mExL + qT * sqrt( diag(W)/(df) * (1/n+1/n) )
  
  mAdL  <- sapply(plastic3[Ad=='L',],mean)
  mAdH  <- sapply(plastic3[Ad=='H',],mean)
  infAd <- mAdH-mAdL - qT * sqrt( diag(W)/(df) * (1/n+1/n) )
  centerAd = mAdH-mAdL
  supAd <- mAdH-mAdL + qT * sqrt( diag(W)/(df) * (1/n+1/n) )
  
  IC2 = list(ExH_ExL=cbind(inf=infEx, center=centerEX, sup=supEx),
             AdH_AdL=cbind(inf=infAd,center=centerAd,sup=supAd))
  IC2
}

# 16)LDA lab8 #####################################################################

## 16.1)Pre-processing ##############################################################

#Fare jittering SE NECESSARIO: 
set.seed(1)
iris2 <- iris2 + cbind(rnorm(150, sd=0.025)) 

indx_a = which(group=='A')
indx_b = which(group=='B')

nA <- length(A)
nB <- length(B)
n <- nA + nB
# Prior probabilities (estimated from the data, no prior knowledge)
PA <- nA / n
PB <- nB / n
### Recall: the classification region is obtained by comparing pA*f.A and pB*f.B
g = 2 #due gruppi
MA <- mean(Infg[A])
MB <- mean(Infg[B])
SA <- var(Infg[A])
SB <- var(Infg[B])
S <- ((nA-1) * SA + (nB-1) * SB) / (nA + nB - g) # pooled estimate
##16.2) Ipotesi #####################################################################
# Assumptions:
# 1) if L=i, X.i ~ N(mu.i, sigma.i^2), i=A,B
# 2) sigma.A=sigma.B
# 3) c(A|B)=c(B|A) (equal misclassification costs)

shapiro.test(cyto[indx_a,1])
var.test(cyto[indx_a,1],cyto[indx_b,1])
bartlett.test(kimono$value, kimono$gruppo)

{#Aggiustare i prior se costi non uguali
# misclassification costs
c.vf <- 10 #costo di giudicare vero un falso
c.fv <- 0.05 #costo di giudicare falso un vero
#prior probabilities
pf <- 0.001 #prior del falso
pt <- 1-0.001 #prior del vero
# Prior modified to account for the misclassification costs
prior.c <- c(pt*c.fv/(c.vf*pf+c.fv*pt), pf*c.vf/(c.vf*pf+c.fv*pt))
prior.c
}

## 16.3)Analisi #####################################################################
{#Call caso univariato (ma la sintassi può essere estesa al multivar!)
cyto.lda <- lda(group ~ Infg) #prior prob estimated from sample
cyto.lda.with.prior <- lda(group ~ Infg, prior=c(0.95,0.05))
cyto.lda
}

{#Call variante con iris2 dataset solo con variabili e species.name vettore factor 
species.name <- factor(Species, labels=c('setosa','versicolor','virginica'))
lda.iris <- lda(iris2, species.name)
}

# In particular:
# - coefficients of linear discriminants: versors of the canonical directions
# [to be read column-wise]
# - proportion of trace: proportion of variance explained by the corresponding
# canonical direction


## 16.4)Prediction ##################################################################

{# Caso univariato: Posterior probability and classification for x new dataframe
x <- data.frame(Infg=seq(-10, 35, 0.5))
# The command predict() returns a list containing (see the help of predict.lda):
# - the class associated with the highest posterior probability
predict(cyto.lda, x)$class
# - the posterior probabilities for the classes
predict(cyto.lda, x)$posterior
# - in lda: the coordinates of the canonical analysis of Fisher
# (Fisher's discriminant scores)
predict(cyto.lda, x)$x
}

{#Caso multivar e calcolo APER con empirical freq
  Lda.iris <- predict(lda.iris, iris2) 
  Lda.iris$class
  table(class.true=species.name, class.assigned=Lda.iris$class)
  errors <- (Lda.iris$class != species.name)
  APER <- sum(errors)/length(species.name)
  APER
  # Remark: this is correct only if we estimate the prior with the empirical
  # frequencies! Otherwise:
  # prior <- c(1/3,1/3,1/3)
  # G <- 3
  # misc <- table(class.true=species.name, class.assigned=Lda.iris$class)
  # APER <- 0
  # for(g in 1:G)
  # APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior[g]
}

{#Estimate of AER by cross validation MANUALLY
  errors_CV <- 0
  prior_probs #da definire
  for(i in 1:n){
    LdaCV.i = lda(iris2[-i,], species.name[-i], prior=prior_probs)
    errors_CV <- errors_CV + as.numeric(predict(LdaCV.i,iris2[i,])$class != species.name[i])
  }
  errors_CV
  AERCV <- sum(errors_CV)/length(species.name)
  AERCV
}

{#Estimate of AER by cross validation WITH R
  LdaCV.iris <- lda(iris2, species.name, CV=TRUE) # specify the argument CV
  table(class.true=species.name, class.assignedCV=LdaCV.iris$class)
  errorsCV <- (LdaCV.iris$class != species.name)
  errorsCV
  AERCV <- sum(errorsCV)/length(species.name)
  AERCV
}

## 16.4) Plot partition ##############################################################

{# Plot the partition induced by LDA
x11()
plot(iris2, main='Iris Sepal', xlab='Sepal.Length', ylab='Sepal.Width', pch=20)
points(iris2[i1,], col='red', pch=20)
points(iris2[i2,], col='green', pch=20)
points(iris2[i3,], col='blue', pch=20)
legend("topright", legend=levels(species.name), fill=c('red','green','blue'), cex=.7)
points(lda.iris$means, pch=4,col=c('red','green','blue') , lwd=2, cex=1.5)
x <- seq(min(iris[,1]), max(iris[,1]), length=200)
y <- seq(min(iris[,2]), max(iris[,2]), length=200)
xy <- expand.grid(Sepal.Length=x, Sepal.Width=y)
z <- predict(lda.iris, xy)$post # these are P_i*f_i(x,y)
z1 <- z[,1] - pmax(z[,2], z[,3]) # P_1*f_1(x,y)-max{P_j*f_j(x,y)}
z2 <- z[,2] - pmax(z[,1], z[,3]) # P_2*f_2(x,y)-max{P_j*f_j(x,y)}
z3 <- z[,3] - pmax(z[,1], z[,2]) # P_3*f_3(x,y)-max{P_j*f_j(x,y)}
17
# Plot the contour line of level (levels=0) of z1, z2, z3:
# P_i*f_i(x,y)-max{P_j*f_j(x,y)}=0 i.e., boundary between R.i and R.j
# where j realizes the max.
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)
contour(x, y, matrix(z3, 200), levels=0, drawlabels=F, add=T)
}

{#Plot in 3d 
  library(rgl)
  library(mvtnorm)
  open3d()
  points3d(iris2[i1,1], iris2[i1,2], 0, col='red', pch=15)
  points3d(iris2[i2,1], iris2[i3,2], 0, col='green', pch=15)
  points3d(iris2[i3,1], iris2[i2,2], 0, col='blue', pch=15)
  surface3d(x,y,matrix(dmvnorm(xy, m1, Sp) / 3, 50), alpha=0.4, color='red')
  surface3d(x,y,matrix(dmvnorm(xy, m2, Sp) / 3, 50), alpha=0.4, color='green', add=T)
  surface3d(x,y,matrix(dmvnorm(xy, m3, Sp) / 3, 50), alpha=0.4, color='blue', add=T)
  box3d()
}

# 17)QDA lab8 #####################################################################

## 17.1) Pre-processing ##############################################################
#Fare jittering SE NECESSARIO: 
set.seed(1)
iris2 <- iris2 + cbind(rnorm(150, sd=0.025)) 

## 17.2)Ipotesi #####################################################################
# Assumptions:
# 1) if L=i, X.i ~ N(mu.i, sigma.i^2), i=A,B
# 2) c(A|B)=c(B|A) (equal misclassification costs)

# Prior modified to account for the misclassification costs
prior.c <- c(pt*c.fv/(c.vf*pf+c.fv*pt), pf*c.vf/(c.vf*pf+c.fv*pt))
prior.c

## 17.3)Analisi #####################################################################
qda.iris <- qda(iris2, species.name)
qda.iris

## 17.4) Prediction ##################################################################
Qda.iris <- predict(qda.iris, iris2)

{#Calcolo APER con prior stimati
  table(class.true=species.name, class.assigned=Qda.iris$class)
  errorsq <- (Qda.iris$class != species.name)
  errorsq
  APERq <- sum(errorsq)/length(species.name)
  APERq
}

{#Calcolo APER con prior NOTI (NON stimati)
prior #da definire
G #numero di gruppi
misc <- table(class.true=species.name, class.assigned=Lda.iris$class)
APER <- 0
for(g in 1:G)
  APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior[g]
}

{#Calcolo AER con cross val USANDO R 
  QdaCV.iris <- qda(iris2, species.name, CV=T)
  table(class.true=species.name, class.assignedCV=QdaCV.iris$class)
  errorsqCV <- (QdaCV.iris$class != species.name)
  errorsqCV
  AERqCV <- sum(errorsqCV)/length(species.name)
  AERqCV
}

## 17.5) Plot partition ##############################################################

{#Caso multiv
  x11()
  plot(iris2, main='Iris Sepal', xlab='Sepal.Length', ylab='Sepal.Width', pch=20)
  points(iris2[i1,], col='red', pch=20)
  points(iris2[i2,], col='green', pch=20)
  points(iris2[i3,], col='blue', pch=20)
  legend("topright", legend=levels(species.name), fill=c('red','green','blue'))
  points(qda.iris$means, col=c('red','green','blue'), pch=4, lwd=2, cex=1.5)
  x <- seq(min(iris[,1]), max(iris[,1]), length=200)
  y <- seq(min(iris[,2]), max(iris[,2]), length=200)
  xy <- expand.grid(Sepal.Length=x, Sepal.Width=y)
  z <- predict(qda.iris, xy)$post
  z1 <- z[,1] - pmax(z[,2], z[,3])
  z2 <- z[,2] - pmax(z[,1], z[,3])
  z3 <- z[,3] - pmax(z[,1], z[,2])
  contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)
  contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)
  contour(x, y, matrix(z3, 200), levels=0, drawlabels=F, add=T)
}

{#Plot 3d 
  open3d()
  points3d(iris2[i1,1], iris2[i1,2], 0, col='red', pch=15)
  points3d(iris2[i2,1], iris2[i3,2], 0, col='green', pch=15)
  points3d(iris2[i3,1], iris2[i2,2], 0, col='blue', pch=15)
  surface3d(x,y,matrix(dmvnorm(xy, m1, S1) / 3, 50), alpha=0.4, color='red')
  surface3d(x,y,matrix(dmvnorm(xy, m2, S2) / 3, 50), alpha=0.4, color='green', add=T)
  surface3d(x,y,matrix(dmvnorm(xy, m3, S3) / 3, 50), alpha=0.4, color='blue', add=T)
  box3d()
  
}

# 18)KNN lab8 #####################################################################

## 18.1) Ipotesi #####################################################################
#Non ha molto senso con spazi di dimensioni grandi -> curse of dimensionality

## 18.2) Pre-processing ##############################################################
set.seed(1)
iris2 <- iris2 + cbind(rnorm(150, sd=0.025)) # jittering

## 18.3) Analisi&Prediction ##########################################################
x <- data.frame(Infg=seq(-10, 35, 0.5))
cyto.knn <- knn(train = Infg, test = x, cl = group, k = 3, prob=T)
attributes(cyto.knn) #mi da levels,class(cioè?) e prob

cyto.knn.class <- (cyto.knn == 'B')+0 
cyto.knn.B <- ifelse(cyto.knn.class==1, 
                     attributes(cyto.knn)$prob, 
                     1 - attributes(cyto.knn)$prob)

{#Scelta di k con plot 
  # let's change k
  x11(width = 28, height = 21)
  par(mfrow=c(3,4))
  for(k in 1:12)
  {
    cyto.knn <- knn(train = Infg, test = x, cl = group, k = k, prob=T)
    cyto.knn.class <- (cyto.knn == 'B')+0 
    cyto.knn.B <- ifelse(cyto.knn.class==1, attributes(cyto.knn)$prob, 1 - attributes(cyto.knn)$prob)
    
    plot(x[,1], cyto.LDA.B, type='l', col='red', lty=2, xlab='x', ylab='estimated posterior', main=k)
    points(x[,1], cyto.knn.B, type='l', col='black', lty=1, lwd=2)
    abline(h = 0.5)
  }
}

## 18.4) Plot partition ##############################################################
{#Plot partizione 2d
k <- 7
x11()
plot(iris2, main='Iris.Sepal', xlab='Sepal.Length', ylab='Sepal.Width', pch=20)
points(iris2[i1,], col=2, pch=20)
points(iris2[i3,], col=4, pch=20)
points(iris2[i2,], col=3, pch=20)
legend("topright", legend=levels(species.name), fill=c(2,3,4))
x <- seq(min(iris[,1]), max(iris[,1]), length=200)
y <- seq(min(iris[,2]), max(iris[,2]), length=200)
xy <- expand.grid(Sepal.Length=x, Sepal.Width=y)
iris.knn <- knn(train = iris2, test = xy, cl = iris$Species, k = k)
z <- as.numeric(iris.knn)
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5), drawlabels=F, add=T)
}


# 19)FISHER DISCRIMINANT ANALYSIS lab8 #################################################
### Let's change viewpoint: we look for the directions that highlight
### the discrimination among groups
### -> we look for the canonical directions

## 19.1) Ipotesi #####################################################################
# Remark. Assumptions: homogeneity of the covariance structure
# [we relax the normality assumption]

## 19.2) Analisi #####################################################################
#Usando il comando della lda  SETTANDO PRIOR UNIFORME,
# -nella parte di "Coefficients of linear discriminants ci sono (leggendo per colonne) le coordinate canoniche
# -nella parte di "Proportion of trace" ci sono le propor di ???? spiegata in quelle direzioni
lda.iris <- lda(iris2, species.name)
lda.iris
a1
a2 #li trascrivo dai numeri che vedo nel output di lda.iris!

## 19.3) Prediction ##################################################################

{### How do I classify a new observation?
x.new <- c(5.85, 2.90)
# compute the canonical coordinates
cc.new <- c(x.new%*%a1, x.new%*%a2)
# compute the distance from the means
dist.m <- c(d1=sqrt(sum((cc.new-cc.m1)^2)),
            d2=sqrt(sum((cc.new-cc.m2)^2)),
            d3=sqrt(sum((cc.new-cc.m3)^2)))
# assign to the nearest mean
which.min(dist.m)
}

## 19.4) Plot prediction #############################################################
{#We plot the partition generated by the canonical coordinates
  color.species <- species.name
  levels(color.species) <- c('red','green','blue')
  
  x11()
  plot(cc1.iris, cc2.iris, main='Fisher discriminant analysis', xlab='first canonical coordinate', ylab='second canonical coordinate', pch=20, col=as.character(color.species))
  legend("topleft", legend=levels(species.name), fill=c('red','green','blue'), cex=.7)
  
  points(cc.m1[1], cc.m1[2], pch=4,col='red' , lwd=2, cex=1.5)
  points(cc.m2[1], cc.m2[2], pch=4,col='green' , lwd=2, cex=1.5)
  points(cc.m3[1], cc.m3[2], pch=4,col='blue' , lwd=2, cex=1.5)
  
  x.cc  <- seq(min(cc1.iris),max(cc1.iris),len=200)
  y.cc  <- seq(min(cc2.iris),max(cc2.iris),len=200)
  xy.cc <- expand.grid(cc1=x.cc, cc2=y.cc)
  
  z  <- cbind( sqrt(rowSums(scale(xy.cc,cc.m1,scale=FALSE)^2)), sqrt(rowSums(scale(xy.cc,cc.m2,scale=FALSE)^2)), sqrt(rowSums(scale(xy.cc,cc.m3,scale=FALSE)^2)))
  z1.cc <- z[,1] - pmin(z[,2], z[,3])    
  z2.cc <- z[,2] - pmin(z[,1], z[,3])    
  z3.cc <- z[,3] - pmin(z[,1], z[,2])
  
  contour(x.cc, y.cc, matrix(z1.cc, 200), levels=0, drawlabels=F, add=T)
  contour(x.cc, y.cc, matrix(z2.cc, 200), levels=0, drawlabels=F, add=T)
  contour(x.cc, y.cc, matrix(z3.cc, 200), levels=0, drawlabels=F, add=T)
  
  dev.off()
}

# 20)SVM lab8 #####################################################################

##20.1) Ipotesi #####################################################################
#Plotta per vedere come sono i dati e guarda se problema è linearmente separabile.
# 1)Se è lin sep, mettendo costo alto impongo hard problem (no overlap)
# 2) Dal plot si capisce quale versione del kernel usare!Varianti possibili sono:
# 'linear', 'radial','polynomial','sigmoid'

## 20.2) Analisi #####################################################################

{#
dat <- data.frame(x=x, y=as.factor (y)) #NB: x must be a matrix , not a dataframe!
svmfit <- svm(y~., data=dat , kernel ='linear', cost =10, scale =FALSE )
#svmfit <- svm(y~., data=dat , kernel ='radial', cost =10, gamma = 1, scale =FALSE ) nel caso radial c'è gamma come parametro aggiuntivo
summary(svmfit)
#From the summary we can read how many support vects are selected in the different groups!
svmfit$index #gives the indexs of the support vectors
}

{#Selezione del parametro di costo con cross val  con kernel='linear'
  
  # To set the parameter C we can use the function tune(),
  # which is based on cross-validation (10-fold)
  
  set.seed (1)
  tune.out <- tune(svm,y~.,data=dat ,kernel = 'linear',
                   ranges =list(cost=c(0.001 , 0.01, 0.1, 1,5,10,100) ))
  summary(tune.out)
  bestmod <- tune.out$best.model
  bestmod$cost #per vedere con quale costo ho ottenuto il modello migliore
  summary(bestmod)
}

{#Scelta costo e gamma con cross val con kernel='radial'
  set.seed (1)
  tune.out <- tune(svm , y~., data=dat[train ,], kernel ='radial',
                   ranges =list(cost=c(0.1 ,1 ,10 ,100 ,1000),
                                gamma=c(0.5,1,2,3,4) ))
  summary(tune.out)
  bestmod <- tune.out$best.model
  bestmod$cost #per vedere con quale costo ho ottenuto il modello migliore
  bestmod$gamma 
  summary(bestmod)
}

## 20.3) Prediction ##################################################################

{#Prediction su un test set 
  testdat <- data.frame(x=xtest , y=as.factor (ytest))
  ypred <- predict(bestmod,testdat)
  table(true.label=testdat$y, assigned.label =ypred )
}

{#Miscassification error (su test o training set)
  
  #-train sono gli indici delle osservazioni usate nel training set
  #-gli indici non in train sono quelli delle osservaz per il test set
  #-"y" è il nome della colonna del gruppo
  
  #Misclassification error on the training set
  table(true=dat[train ,"y"], pred=predict (svmfit ,
                                            newdata =dat[train ,]))
  # Misclassification error on the test set
  table(true=dat[-train ,"y"], pred=predict (svmfit ,
                                             newdata =dat[-train ,]))
  
  # Miscl_error = Dalla tabella, somma le misclass (ossia fuori diag) e dividi per la somma di
  #               tutti gli elementi della tabella 
}

## 20.4)Plot prediction #############################################################

{#Plot di default (ma ha gli assi scambiati!)
  x11()
  par(mfrow=c(1,2))
  plot(svmfit , dat, col =c('salmon', 'light blue'), pch=19, asp=1)
}

{#Plot con griglia colorata di classificazione
  n.g <- 100
  xgrid <- expand.grid(x.1=seq(from=range(dat$x.1)[1],to=range(dat$x.1)[2],length=n.g),
                       x.2=seq(from=range(dat$x.2)[1],to=range(dat$x.2)[2],length=n.g))
  ygrid <- predict(svmfit,xgrid)
  x11()
  plot(xgrid,col=c("red","blue")[as.numeric(ygrid)],pch=20,cex=.2)
  points(x,col=c("red","blue")[as.numeric(dat$y)],pch=19)
  points(x[svmfit$index,],pch=5,cex=2)
}

{#Plot con linea di separazione
  x11()
  plot(x,col=c("red","blue")[as.numeric(dat$y)],pch=19)
  contour(seq(from=range(dat$x.1)[1],to=range(dat$x.1)[2],length=n.g),
          seq(from=range(dat$x.2)[1],to=range(dat$x.2)[2],length=n.g),
          matrix(as.numeric(ygrid),n.g,n.g),level=1.5,add=TRUE,
          drawlabels=F)
}


# 21) HIERARCHICAL CLUSTERING lab9 #############################################
## 21.1) Pre-processing ########################################################

#We can view the data through the image command and see if they seem to be ordered
# If yes, we can destruct the order resampling :

x11()
par(mfrow=c(1,3))
image(1:n,1:n,as.matrix(dist.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j' )
image(1:n,1:n,as.matrix(dist.c), main='metrics: Canberra', asp=1, xlab='i', ylab='j' )
image(1:n,1:n,as.matrix(dist.m), main='metrics: Manhattan', asp=1, xlab='i', ylab='j' )

disordered_indxs <- sample(n) #estraggo i 150 elementi in ordine casuale
iris4 <- iris4[disordered_indxs,]

## 21.2) Analisi ###############################################################

d = dist(dataset, method='euclidean') #'manhattan','canberra'
clusters = hclust(d, method='single') #'average','complete', 'ward'

clusters$merge # ordina per ordine di aggregazione
clusters$height # distanza a cui si ha aggregazione
clusters$order # ordering that allows to avoid intersections in the dendrogram

plot(clusters, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='') #plotta il dendogram
rect.hclust(clusters, k) #plotta sul dendogram i rettangoli dei k cluster

clusters_tagliati = cutree(clusters, k=2) #taglia il dendogram quando arriva a k cluster

plot(dataset, col=ifelse(cluster_tagliati==1,'red','blue'), pch=19)

#Se si hanno le etichette, si possono contare gli errori
table(label.true = species.name[disordered_indxs], label.cluster = clusters_tagliati) #NB: si usa species.name[disordered_indxs] se si è fatto resampling all'inizio
table(clusters_tagliati) #mi fa vedere la dimensione dei vari clusters

coph_matr = cophenetic(clusters) #usare su cluster NON tagliati

coph_coef = cor(dist, coph_matr) #calcola il cophenetic coeff 

#Confronto matrice delle distanze prima e dopo cluster
image(as.matrix(dx), asp=1, main='not ordered x' )
dx <- dist(x)
hcx<- hclust(dx, method='single')
image(as.matrix(dx)[hcx$order,hcx$order], asp=1, main='reordered x' )


# 22) K-MEANS ##################################################################

# end) LIBRARY ######################################################################

library(plotrix) #usata nel lab1 per grafici torta3d
library(rgl) # lab1 per grafici3d e vari
library(mvtnorm) #lab3 per generare dati
library(car) #lab3 per plot
library(MASS) #per lda
library(class) #per knn
library(e1071) #per SVM
library (ISLR) #lab8 alla fine,non so per cosa

