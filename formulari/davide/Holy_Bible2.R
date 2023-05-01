#### Linear models ###################################
### Model:
### distance = beta_0 + beta_1 * speed + beta_2 * speed^2 + Eps
### Assumptions:
## 1) Parameter estimation: E(Eps) = 0  and  Var(Eps) = sigma^2 
## 2) Inference:            Eps ~ N(0, sigma^2)
##### 1) Estimate of the parameters
##### Assumptions: E(Eps) = 0  and  Var(Eps) = sigma^2 
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

# funzione di day
fun=1+cos((4*pi*day)/365)

############################# Inference on the parameters #####
### Assumption: Eps ~ N(0, sigma^2)
### Test (Fisher):
shapiro.test(residuals(fm))
# H0: (beta1, beta2) == (0, 0) vs H1: (beta1, beta2) != (0, 0)
linearHypothesis(fm, rbind(c(0,1,0), c(0,0,1)), c(0,0)) 
# i vettori sono i coeff che metto prima dei parametri: tipo sopra sto 
# confrontando b1 (primo vettore) con b2 (secondo)
summary(fm)
# i test li puoi fare a mano usando la t statistics (solo Oat!)
# linearHypothesis serve per testarli both contemporaneamente

p <- 2  # number of tested coefficients
r <- 2  # number of regressors

# Confidence region: PER I PARAMETRI CIOE I BETA NON PER LE OBS
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

############################# Confidence intervals for the mean & prediction (new obs) ####
##### Assumption: Eps ~ N(0, sigma^2)
# Command predict()
Z0.new <- data.frame(speed1=10, speed2=10^2) # questo è il nuovo dato
# Conf. int. for the mean
Conf <- predict(fm, Z0.new, interval='confidence', level=1-0.05)  
Conf  # intervallo sulla media campionaria delle nuove oobs
# Pred. int. for a new obs
Pred <- predict(fm, Z0.new, interval='prediction', level=1-0.05)  
Pred  # intervallo sulle nuove obs, piu grande
# Questo mi da la stima li
Resp <- predict(fm, Z0.new, type = 'response') 

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

#in alternativa usa matplot: dopo x11 e plot
matplot(Z0,Conf,add=T,type='l',col=c('black','blue','blue'),lwd=2,lty=2)
matplot(Z0,Pred,add=T,type='l',col=c('black','green','green'),lwd=2,lty=2)
legend('topleft', legend=c('regression line','confidence intervals','prediction intervals'),
       col=c('black','blue','green'), lwd=2, cex=0.85)

# Bonferroni intervals sulla media e sulla var
k <- 2 # numero di intervalli
alpha <- .1
n <- dim(index)[1]
r <- 3 # numero di beta
Z0   <- data.frame(D=1, Anno=1800)
ICBmean <- predict(fitB, Z0, interval='confidence',level=1-alpha/k) 
e <- residuals(fitB)
ICBvar <- data.frame(L=t(e)%*%e/qchisq(1-alpha/(2*k),n-(r+1)),
                     U=t(e)%*%e/qchisq(alpha/(2*k),n-(r+1)))

# valore che massimizza la predizione
tm <- which.max(fm2$coefficients[2]*1+cos((4*pi*day)/365) + fm2$coefficients[3]*d_wend)
t.max <- data.frame(fun=1+cos((4*pi*tm)/365), d_wend=d_wend[tm])
Pred <- predict(fm2, t.max, interval='prediction', level=1-0.05)  
Pred
# nei modelli polinomiali/continui fai test sulla derivata sul valore fornito (61)
A <- c(0,1,2*61)
b <- 0
linearHypothesis(fm2,A,b)

############################# Verify assumptions ####
x11()
par(mfrow=c(2,2))
plot(fm)  # fa tutti i plot utili
shapiro.test(residuals(fm))
# If the pvalue of the Shapiro-Wilk Test is greater than 0.05, the data is normal

# se non hai normalita: prova trasformazione o togli outliers
which(fm$residuals>1.5) # 1 16 24
data=cbind(temp1,temp2,lun)
data=data[-which(fm$residuals>1.5),]

#plot della retta di regressione nel caso 2D
result <- lm(Y ~ X)
coef <- result$coef
plot(X, Y)
abline(coef[1],coef[2], lwd=2,col='red')

### Example 4: Earthquakes
# Model:
# depth = beta0 + beta1*lat + beta2*long + eps
fit  <- lm(depth ~ lat + long, data=Q) # modo per chiamare direttamente dal dataset Q

############################# dummy variable ####
# creare velocemente dummy (2 factor)
d_wend <- rep(0,length(day.of.the.week))
d_wend[which(day.of.the.week=='weekend')] <- 1
# per 3 factor usa 2 dummy: ripeti 2 volte la procedura sopra

# stima i coeff del modello (con dummy) (teoria spiegata sotto)
fm <- lm(waiting.time ~ fun*d_wend) 
summary(fm)
cf=coefficients(fm)
coeff=c(cf[1],cf[1]+cf[3],cf[2],cf[2]+cf[4],sqrt(sum(residuals(fm)^2)/fm$df))
coeff

# dummy coi cluster
dummy <- clusterw - 1   # 0 = red, clusterw arriva da cutree, 1 = green
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


#### Problem of collinearity ##############################
# solo push the button
# criprendo cars:
n          <- dim(cars)[[1]]
distance   <- cars$dist
speed1     <- cars$speed
speed2     <- cars$speed^2
fm <- lm(distance ~ speed1 + speed2)  #x+x^2
summary(fm) 
vif(fm) #alto, significa che ho collinearity: sopra 5 o 10 la var è sus

############################# PCA ##############################################
speed.pc <- princomp(cbind(speed1,speed2), scores=TRUE)
summary(speed.pc)
speed.pc$load
sp1.pc <- speed.pc$scores[,1]
sp2.pc <- speed.pc$scores[,2]
# Now we estimate the model by inserting the PCs instead of the original X
fm.pc <- lm(distance ~ sp1.pc + sp2.pc)
summary(fm.pc) 
# vedo che solo la prima principal component è significativa (pvalue basso)
# => dimensionality reduction: fm.pc <- lm(distance ~ sp1.pc)

############################# ridge regression ##################################
lambda <- .5
fit.ridge <- lm.ridge(distance ~ speed1 + speed2, lambda = lambda)
# Note: to fit the model, R automatically centers X and Y 
# with respect to their mean.
coef.ridge <- coef(fit.ridge)
yhat.lm <- cbind(rep(1,n), speed1, speed2)%*%coef(fm)  # LM fitted values
yhat.r <- cbind(rep(1,n), speed1, speed2)%*%coef.ridge # ridge fitted values

# cerco lambda opt
lambda.c <- seq(0,10,0.01)
fit.ridge <- lm.ridge(distance ~ speed1 + speed2, lambda = lambda.c)
lambda.opt <- lambda.c[which.min(fit.ridge$GCV)]
coef.ridge <- coef(fit.ridge)[which.min(fit.ridge$GCV),]

############################# lasso ###########################################
# Build the matrix of predictors
x <- model.matrix(distance~speed1+speed2)[,-1] #model matrix utile con le cat
# Build the vector of response
y <- distance  # devo per forza passargli x e y cosi, un po scomoda

# Let's set a grid of candidate lambda's for the estimate
lambda.grid <- 10^seq(5,-3,length=100) # piu utile prendere punti 
# equidistribuiti in scala logaritmica
fit.lasso <- glmnet(x,y, lambda = lambda.grid) # default: alpha=1 -> lasso 
# [note: if alpha=0 -> ridge regression]

# Let's set lambda via cross validation
cv.lasso <- cv.glmnet(x,y,alpha=1,nfolds=10,lambda=lambda.grid) # default: 10-fold CV
bestlam.lasso <- cv.lasso$lambda.min

plot(cv.lasso)
abline(v=log(bestlam.lasso), lty=1)

# Get the coefficients for the optimal lambda
coef.lasso <- predict(fit.lasso, s=bestlam.lasso, type = 'coefficients')

############################# Multiple linear regression ##########################
data <- read.table('concrete.txt', header=T) #esempio pratico

### Ridge and Lasso regression with glmnet (strada alternativa)

x <- model.matrix(Y ~ X1 + X2 + X3 + X4)[,-1] # matrix of predictors
y <- Y # vector of response
lambda.grid <- 10^seq(5,-3,length=50)

# Ridge regression
fit.ridge <- glmnet(x,y, lambda = lambda.grid, alpha=0) # alpha=0 -> ridge 

x11()
plot(fit.ridge,xvar='lambda',label=TRUE, col = rainbow(dim(x)[2]))
legend('topright', dimnames(x)[[2]], col =  rainbow(dim(x)[2]), lty=1, cex=1)
abline(v=log(bestlam.ridge)) #aggiungi questa linea solo dopo che trovi bestmod

# Let's set lambda via CV
set.seed(1)
cv.ridge <- cv.glmnet(x,y,alpha=0,nfolds=3,lambda=lambda.grid)
bestlam.ridge <- cv.ridge$lambda.min

plot(cv.ridge)
abline(v=log(bestlam.ridge), lty=1)
#linea scura: optimal lambda scelto dal crterio
#linea trattegiata: lambda + grande st sono dentro 1 deviazione standard dal optimal lambda
coef.ridge <- predict(fit.ridge, s=bestlam.ridge, type = 'coefficients')

### Lasso regression
fit.lasso <- glmnet(x,y, lambda = lambda.grid, alpha=1) # alpha=1 -> lasso 
#puoi riciclare tutto il codice di prima same plot di prima

# Compare coefficients estimates for LS, Ridge and Lasso
# l2 norm
sqrt(sum((coef(lm(y~x))[-1])^2)) # LS
sqrt(sum((coef.ridge[-1])^2))    # ridge
# l1 norm
sum(abs(coef(lm(y~x))[-1])) # LS
sum(abs(coef.lasso[-1]))    # lasso
# l2 norm vince ridge l1 vicne lasso

### Variable selection
linearHypothesis(result, rbind(c(0,0,1,0,0),c(0,0,0,1,0)),c(0,0))
# modo per togliere 2 covariate in un colpo solo

# Best Subset Selection (exhaustive search)
regfit.full <- regsubsets(Salary~., data=Hitters, nvmax=19) #num max cov
reg.summary <- summary(regfit.full)
reg.summary$which
reg.summary$rsq   # r-squared
reg.summary$adjr2 # adjusted r-squared
reg.summary$rss   # residual sum of squares

x11(height=7,width=14)
par(mfrow=c(1,3))
plot(reg.summary$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
plot(reg.summary$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
plot(reg.summary$rss,xlab="Number of Variables",ylab="RSS",type="b")

which.max(reg.summary$adjr2)
coef(regfit.full,11) #estrae i coeff

# Forward and Backward Stepwise Selection
regfit.fwd <- regsubsets(Salary~.,data=Hitters,nvmax=19,method="forward")
# adesso aggiungo e basta, le covariate una volta inserite non le tolgo piu 
summary(regfit.fwd)

#al contrario: perto con tutte e le elimino oat
regfit.bwd <- regsubsets(Salary~.,data=Hitters,nvmax=19,method="backward")
summary(regfit.bwd)
#stessi grafici di prima


#### Linear mixed models ####################################
library(nlmeU) ## --> for the dataset
library(nlme)  ## --> for models implementation
library(corrplot)
library(lattice)
library(plot.matrix)
library(ggplot2)
library(insight)
library(lattice)
library(lme4)
#dataset armd
## sample means across time and treatment
flst <- list(armd$time.f, armd$treat.f)
tMn <- tapply(armd$visual, flst, FUN = mean)
tMn

# 1. Linear Models with homogeneous and independent errors: we start by considering all 
#    observations as independent, with homogeneous variance

# LM: VISUAL_it = b_0t + b1 × VISUAL_0i + b_2t × TREAT_i + e_it

lm1.form <- lm(visual ~ -1 + visual0 + time.f + treat.f:time.f, data = armd )
# con il -1 rimuovo l'intercetta del modello
summary(lm1.form)
sig=sum(residuals(lm1.form)^2)/lm1.form$df

# variance-covariance matrix of Y  --> it is a diagonal matrix with a value of sig
sig=sum(residuals(lm1.form)^2)/lm1.form$df
x11()
par(mar = c(4,4,4,4))
plot(diag(x=sig,nrow=30, ncol=30), main='Variance-covariance matrix of Y')

# vedo come cambiano i residui per ogni paziente
colori =rainbow(length(unique(armd$subject)))
boxplot(lm1.form$residuals ~ armd$subject, col=colori,
        xlab='Subjects', ylab='Residuals', main ='Distribution of residuals across patients') 
# no omoschedasticity: non tutti i boxplot hanno media 0, e hanno varianze diverse

############################# Linear models with heteroscedastic and independent errors ####
fm9.1 <- gls(visual ~ -1 + visual0 + time.f + treat.f:time.f,  # the same as before
             weights = varIdent(form = ~1|time.f), # Var. function; <delta, stratum>-group
             data = armd)  #la funzione gls inlcude la dipendenza e eterosched
# la var dipende sia dal soggetto che dal tempo: in  questo caso solo etero
summary(fm9.1)$par

fm9.1$modelStruct$varStruct[1,] #Parameter estimates
intervals(fm9.1, which = "var-cov")  ## 95% CI sulla varaiance function
#Parameter estimates:
#    4wks    12wks    24wks    52wks 
#1.000000 1.397600 1.664321 1.880852

fm9.1$sigma #Residual standard error
par(mar = c(4,4,4,4))
plot(diag(x=c(1.000000^2*fm9.1$sigma^2, 1.397600^2*fm9.1$sigma^2, 1.664321^2*fm9.1$sigma^2, 1.880852^2*fm9.1$sigma^2), nrow=30, ncol=30),
     main='Variance-covariance matrix of Y - VarIdent()')

## To formally test the hypothesis that the variances are timepoint specific, 
## we apply the anova() function. The LR test tests the null hypothesis of homoscedasticity.

anova(fm9.1, lm1.form)  ## lm1.form C fm9.1
# p value basso => signifcativo, dunque la media varia tra i gruppi, il treat ha effetto

## 2.2 Option 2: VarPower()
fm9.2 <- update(fm9.1, weights = varPower(form = ~time)) # Var. function; <delta, v_it>-group
summary(fm9.2)  #modello + facile: 4 sigma diversi, uno per ogni time istant come prima
# ma stavolta per stimarli uso solo un parametro, cioe power

fm9.2$modelStruct$varStruct  #Parameter estimates: power=0.2519332 
intervals(fm9.2, which = "var-cov")   ## 95% CI sulla varaiance function

fm9.2$sigma #Residual standard error
par(mar = c(4,4,4,4))
plot(diag(x=c(4^(2*0.2519332)*fm9.2$sigma^2, 12^(2*0.2519332)*fm9.2$sigma^2, 24^(2*0.2519332)*fm9.2$sigma^2, 52^(2*0.2519332)*fm9.2$sigma^2), nrow=30, ncol=30),
     main='Variance-covariance matrix of Y - VarPower')
#4,12,24,52 sono i 4 valori assunti da time

# Test of the variance structure: power of time vs. timepoint-specific variances
anova(fm9.2, fm9.1)  # p-value alto => uguali
AIC(fm9.2, fm9.1)  # --> fm9.2 is better in terms of AIC and parsimony!

## Residual analysis: time.f categorica
plot(fm9.2, resid(., type = "response") ~ fitted(.)) # Raw vs. fitted
bwplot(resid(fm9.2) ~ time.f, pch = "|", data = armd)
# The boxand-whiskers plots clearly show an increasing variance of the residuals.
## Pearson residuals 
plot(fm9.2, resid(., type = "pearson" ) ~ fitted(.)) # Pearson vs. fitted
bwplot( resid(fm9.2, type = "pearson") ~ time.f, # Pearson vs. time.f
        pch = "|", data = armd)  ## this plot illustrate the effect of scaling:
# the variance of the residuals is virtually constant.

############################# Linear models with heteroscedastic and dependent errors ####
## We now modify the model, so that the visual acuity measurements,
## obtained for the same individual, are allowed to be correlated.
## Variogram per time difference 
Vg1 <- Variogram(fm9.2, form = ~ time | subject)
Vg1
plot(Vg1, smooth = FALSE, xlab = "Time difference",ylim=c(0,0.7))

## From this plot we see that correlation decreases with time lag/difference
## A more appropriate structure might be, e.g., an autoregressive process of order 1 AR(1).

## 3.1 Correlation 1: CorCompSym()
lm1.form <- formula(visual ~ -1 + visual0 + time.f + treat.f:time.f )
fm12.1 <- gls(lm1.form, weights = varPower(form = ~time),
              correlation = corCompSymm(form = ~1|subject),
              data = armd)
summary(fm12.1)
intervals(fm12.1, which = "var-cov")
# With the estimates of rho, sigma and delta we can estimate the var-cov matrix
# The marginal variance-covariance structure
fm12.1vcov <- getVarCov(fm12.1, individual = "2")  #estimate of R_i, e.g. i=2
# penso che se cambi individual non cambia un cazzo 

## on the diagonal we have (sig^2)*TIME^(2*power)
## out of the diagonal we have (sig^2)*TIME_1^(power)*TIME_2^(power)*rho
print(cov2cor(fm12.1vcov), corr = TRUE, stdevs = FALSE)  ## Estimate of C_i (correlation matrix)

## Visualization of the marginal variance-covariance matrix of Y
# devi copiare brutalmente fm12.1vcov
R_i =rbind(c( 73.531 , 56.077 , 67.143 , 82.081),
           c(56.077, 130.140 , 89.323, 109.200),
           c(67.143,  89.323, 186.560, 130.740),
           c(82.081, 109.200, 130.740, 278.810))

R = matrix(0, nrow=28, ncol=28)
for(i in 0:6){
  R[(i*4+1):(i*4+4),(i*4+1):(i*4+4)] = R_i
}
plot(R)

## Test of independence vs. compound-symmetry correlation structure
anova(fm9.2, fm12.1) # M9.2 C M12.1, il test ha effeto

## 3.2 Correlation 2: AR(1)
fm12.2 <- update(fm9.2, 
                 correlation = corAR1(form = ~tp|subject),
                 data = armd)
# poi stesse sbaraccate di prima compreso anova
## 3.3 Correlation 3: general correlation structure, SALTATO DA QUI IN POI
fm12.3 <- update(fm12.2, correlation = corSymm(form = ~tp|subject),  ## the variance function is still VarPower()
                 data = armd)

## Model-Fit Diagnostics
plot(fm12.3)  # residuals, if the plot reveals a pattern:
plot(fm12.3,  # residuals for each time istant
     resid(., type = "p") ~ fitted(.) | time.f)

############################# Random intercept, homoscedastic residuals ####################
fm16.1mer <- lmer(visual ~ visual0 + time * treat.f + (1|subject), #(1|subject) 
                  data = armd)  #è la random intercept at a subject level
#  Groups   Name        Variance Std.Dev.
# subject  (Intercept) 80.61    8.978    lei è d11
# Residual             74.43    8.628    lei la classifica sigmasquared
summary(fm16.1mer)
confint(fm16.1mer,oldNames=TRUE)

## Var-Cov matrix of fixed-effects
vcovb <- vcov(fm16.1mer) #cov
corb <- cov2cor(vcovb) #corr

## Var-Cov matrix of random-effects and errors
print(vc <- VarCorr(fm16.1mer), comp = c("Variance", "Std.Dev."))
# le estraggo:
sigma2_eps <- as.numeric(get_variance_residual(fm16.1mer))
sigma2_b <- as.numeric(get_variance_random(fm16.1mer))

## Let's compute the conditional and marginal var-cov matrix of Y
sgma <- summary(fm16.1mer)$sigma
A <- getME(fm16.1mer, "A") # A  --> N x n, A represents the D (not italic)
I.n <- Diagonal(ncol(A)) # IN  --> n x n
## the conditional variance-covariance matrix of Y (diagonal matrix)
SigmaErr = sgma^2 * (I.n)
SigmaErr[3:6, 3:6]  ## visualization of individual 2

## we visualize the first 20 rows/columns of the matrix
plot(as.matrix(SigmaErr[1:20,1:20]), main = 'Conditional estimated Var-Cov matrix of Y')

## the marginal variance-covariance matrix of Y (block-diagonal matrix)
V <- sgma^2 * (I.n + crossprod(A)) # V = s^2*(I_N+A*A) --> s^2*(I_N) is the error part, s^2*(A*A) is the random effect part
V[3:6, 3:6]  #-> V is a block-diagional matrix, the marginal var-cov matrix

# visualization of the first 20 rows/columns
plot(as.matrix(V[1:20,1:20]), main = 'Marginal estimated Var-Cov matrix of Y')

# Another way to interpret the variance output is to note percentage of the subject variance out 
# of the total, i.e. the Percentage of Variance explained by the Random Effect (PVRE).
# This is also called the intraclass correlation (ICC), because it is also an estimate of the within 
# cluster correlation.

PVRE <- sigma2_b/(sigma2_b+sigma2_eps)
PVRE # 51% is very high!

## visualization of the random intercepts with their 95% confidence intervals
# Random effects: b_0i for i=1,...,234
dotplot(ranef(fm16.1mer, condVar=T))

# The dotplot shows the point and interval estimates for the random effects, 
# ordering them and highlighting which are significantly different from the mean (0)

# Let's now examine standard predictions vs. subject-specific predictions.
# Prediction from mixed model on a test observation from a subject present in the training set:
test.data= data.frame(subject= '234', treat.f='Active', visual0= 63, time = 12)
# 1) Without random effects ->  re.form=NA
predict_no_re <- predict(fm16.1mer, newdata = test.data, re.form=NA)
# 2) With random effects
predict_re <- predict(fm16.1mer, newdata = test.data)

# Prediction from mixed model on a test observation from a subject NOT present in the training set:
test.data= data.frame(subject= '400', treat.f='Active', visual0= 63, time = 12)
# 1) Without random effects ->  re.form=NA
predict_no_re <- predict(fm16.1mer, newdata = test.data, re.form=NA)
# 2) With random effects
predict_re <- predict(fm16.1mer, newdata=test.data, allow.new.levels = T)

# Diagnostic plots 
# 1) Assessing Assumption on the within-group errors
x11()
plot(fm16.1mer)  ## Pearson and raw residuals are the same now

x11()
qqnorm(resid(fm16.1mer))
qqline(resid(fm16.1mer), col='red', lwd=2)

# 2) Assessing Assumption on the Random Effects
x11()
qqnorm(unlist(ranef(fm16.1mer)$subject), main='Normal Q-Q Plot - Random Effects on Intercept')
qqline(unlist(ranef(fm16.1mer)$subject), col='red', lwd=2)

############################# random intercept + slope and homoscedastic residuals ####
fm16.2mer <- lmer(visual ~ visual0 + time * treat.f + (1+time|subject),
                  data = armd, control=lmerControl(optimizer="bobyqa",
                                                   optCtrl=list(maxfun=2e5)))
# slope sia sul tempo che sul soggetto, stesso codice di prima

# PVRE
# In this case the variance of random sigma2_R effects represents the mean random 
# effect variance of the model and is given by
# sigma2_b = Var(b0,b1) = sigma2_b0 + 2Cov(b0,b1)*mean(w) + sigma2_b1*mean(w^2)
# per il resto uguale

# Diagnostic plots , tutto uguale tranne
# 2) Assessing Assumption on the Random Effects
x11()
qqnorm(unlist(ranef(fm16.2mer)$subject[,1]), main='Normal Q-Q Plot - Random Effects on Intercept')
qqline(unlist(ranef(fm16.2mer)$subject[,1]), col='red', lwd=2)

x11()
qqnorm(unlist(ranef(fm16.2mer)$subject[,2]), main='Normal Q-Q Plot - Random Effects on Slope')
qqline(unlist(ranef(fm16.2mer)$subject[,2]), col='red', lwd=2)

############################# diagonal D ######################################
fm16.2dmer <- lmer(visual ~ visual0 + time * treat.f + (1|subject) + (0 + time|subject),
                   data = armd, control=lmerControl(optimizer="bobyqa",
                                                    optCtrl=list(maxfun=2e5)))
# tutto uguale
# PVRE

# In this case the variance of random sigma2_R effects represents the mean random 
# effect variance of the model and is given by
# sigma2_b = Var(b0,b1) = sigma2_b0 + 0 + sigma2_b1*mean(z^2)

############################# LINEAR MIXED MODELS WITH HETEROSCEDASTIC RESIDUALS ####
## fixed-effects formula
lm2.form <- formula(visual ~ visual0 + time + treat.f + treat.f:time ) 
# LMM with homoscedastic residuals
fm16.1 <- lme(lm2.form, random = ~1|subject, data = armd)
# update fm16.1 including heteroscedastic residuals
fm16.2 <- update(fm16.1,
                 weights = varPower(form = ~ time), 
                 data = armd)
summary(fm16.2)
VarCorr(fm16.2) 

## var-cov matrix of the errors (i.e. of Y, conditional to the random effects), that are independent but heteroscedastic 
fm16.2ccov = getVarCov(fm16.2, type = "conditional",  individual = "2")
fm16.2ccov  # diagonale
plot(as.matrix(fm16.2ccov[[1]]), main = expression(paste('Conditional estimated Var-Cov matrix of ', Y[2])))

## var-cov matrix of Y_i
fm16.2cov = getVarCov(fm16.2, type = "marginal", individual = "2")
fm16.2cov # (90.479 = 31.103 + 59.37555; 121.440 = 62.062 + 59.37555; ...)
plot(as.matrix(fm16.2cov[[1]]), main = expression(paste('Marginal estimated Var-Cov matrix of ', Y[2])))

## correlation matrix of Y_i
cov2cor(fm16.2cov[[1]])

## ANALYSIS OF RESIDUALS, niente di nuovo
plot(fm16.2) # Default residual plot of conditional Pearson residuals
# Plots (and boxplots) of Pearson residuals per time and treatment
plot(fm16.2, resid(., type = "pearson") ~ time | treat.f,
     id = 0.05)
bwplot(resid(fm16.2, type = "p") ~ time.f | treat.f, 
       panel = panel.bwplot, # User-defined panel (not shown)
       data = armd)
# Normal Q-Q plots of Pearson residuals 
qqnorm(fm16.2, ~resid(.) | time.f) 

## ANALYSIS OF RANDOM EFFECTS
# Normal Q-Q plots of predicted random effects
qqnorm(fm16.2, ~ranef(.))  
## Computing predictions comparing population average predictions with patient-specific predictions
aug.Pred <- augPred(fm16.2,
                    primary = ~time, # Primary covariate
                    level = 0:1, # fixed/marginal (0) and subj.-spec.(1)
                    length.out = 2) # evaluated in two time instants (4 e 52 wks)

plot(aug.Pred, layout = c(4, 4, 1))
# altri modelli omessi secondo me non fatti

############################# aggiunte di fine capitolo a cazzo #####################
# Random effects: b_0i
ranef(lmm1)

# The dotplot shows the point and interval estimates for the random effects, 
# ordering them and highlighting which are significantly different from the mean (0)

x11()
dotplot(ranef(lmm1))

# in generale in 3.LMM 2.R esercizi utili su lmer
# cioe modelli con Random Intercept & Slope e Linear Mixed Effects Models


#### Spacial Statistics ###########
library(sp)           ## Data management
library(lattice)      ## Data management
library(gstat)        ## Geostatistics
## Define the sample coordinates
coordinates(meuse) <- c('x','y')  #dico a r che ho un dataset speciale
# legge le prime 2 colonne come coordinate, x e y sono i nomi che do alle colonne
# adesso meuse è un dataset speciale, posso fare special plots

# FAI SEMPRE BUBBLE E SCATTER, soprattutto scatter

# bubble plot(obj,zcol,...)
# key.space=location of the key
bubble(meuse,'zinc',do.log=TRUE,key.space='bottom')  # zinc è il nome della var
# le bolle sono le posizioni dei dati, + grande le palle + grande il valore della var

# river meuse
data(meuse.riv)
meuse.lst <- list(Polygons(list(Polygon(meuse.riv)), "meuse.riv"))
meuse.sr <- SpatialPolygons(meuse.lst) # il fiume

# grid for prediction
data(meuse.grid)  # è l'intorno del fiume
is(meuse.grid)
coordinates(meuse.grid) <- c('x','y')
meuse.grid <- as(meuse.grid, 'SpatialPixelsDataFrame')

# plot all together
image(meuse.grid, col = "lightgrey")
plot(meuse.sr, col = "grey", add = TRUE)
plot(meuse, add = TRUE) # x per ogni dato
title('meuse river geostatistical data')

# scatterplot of log(zinc) with respect to distance from the river 
xyplot(log(zinc) ~ sqrt(dist), as.data.frame(meuse)) # vedo la corr col dato spaziale

############################# Variogram Analysis ####
# list of parametric isotropic variogram models
vgm()
# STESSI COMANDI PER UNISOTROPY, SOLO NON METTERE ~ 1 MA ~ D
# sample variogram (binned estimator)
svgm <- variogram(log(zinc) ~ 1, meuse) #solo l'intercetta cioe 1
plot(svgm, main = 'Sample Variogram',pch=19) # isotropy means homogeneity in all directions
# modello isotropico, assumo spatially constant mean,
plot(variogram(log(zinc) ~ 1, meuse, alpha = c(0, 45, 90, 135)),pch=19)
# variogram in piu dimensioni: poi ci sarebbero altri parametri grafici

# dal plot del variogram cerca di capire un po a occhio il modello
# Recall: both spherical and exponential model have a linear behavior near the
#         origin but exponential model has a faster growth than the spherical one

# plot of the final fit
v <- variogram(log(zinc) ~ 1, meuse)
v.fit <- fit.variogram(v, vgm(1, "Sph", 800, 1)) #800 è il range,
# cioe il valore per cui il variogram raggiunge il sill (asymptotic behaviour)
# il primo parametro è il sill, asintoto orizzontale, il quarto è il nuggget
plot(v, v.fit, pch = 19)

# 3 modi per fittare il modello:
attr(v.fit, 'SSErr') # non ho capito che fa
fit.variogram(v, vgm(1, "Sph", 800, 0.06), fit.sills = c(FALSE, TRUE))
# con c(FALSE, TRUE) fitto solo il sill legato allo sph model, non il nugget
# the range parameters can be fixed using argument fit.ranges
fit.variogram.reml(log(zinc)~1, meuse, model=vgm(0.6, "Sph", 800, 0.06))
v.fit # fitto con la likehood

## modeling anisotropy*, modello un var per ogni dir
v.dir <- variogram(log(zinc)~1,meuse,alpha=(0:3)*45) 
v.anis <- vgm(.6, "Sph", 1600, .05, anis=c(45, 0.3))
print(plot(v.dir, v.anis, pch=19))

############################# SPATIAL PREDICTION & KRIGING #############################
## Prediction at a single new location, costruisco nuovo dato
s0.new=data.frame(x=179180, y=330100) # UTM coordinates
coordinates(s0.new)=c('x','y')

plot(s0.new, add = TRUE, col='red', lwd = 2) #aggiungilo a plot all together

# Create a gstat object setting a spherical (residual) variogram
# gstat(g.obj, id, formula, data, model, set,...)
g.tr <- gstat(formula = log(zinc) ~ 1, data = meuse, model = v.fit) # model=variogram

## ordinary kriging
predict(g.tr, s0.new) # questa è la mia point prediction
# variance > 0 (as expected), perche il dato è nuovo
# Estimate the mean:use the argument 'BLUE'
predict(g.tr, s0.new, BLUE = TRUE) #mi da il blue for the mean
# this gives the estimate of the mean under gls, nei modelli lin coincidono 
# stima che non tiene conto della spatial dependence
predict(g.tr, meuse[1,]) # predicta esattamente log di zinc, var=0
# this gives the estimate of the mean (drift component) under gls
predict(g.tr, meuse[1,], BLUE = TRUE) # lui non è esatto! ma per ogni punto ho lam stessa varianza
# prediction over the entire grid
lz.ok <- predict(g.tr, meuse.grid, BLUE = FALSE)
spplot(lz.ok) # noto uno special pattern

#### Universal Kriging
# the hypothesis of spatially constant mean may be too restrictive!

# Create a gstat object setting a spherical (residual) variogram
# gstat(g.obj, id, formula, data, model, set,...)
meuse.gstat <- gstat(id = 'zinc', formula = log(zinc) ~ sqrt(dist),
                     data = meuse, nmax = 50, model=v.fit, set = list(gls=1))
meuse.gstat  # formula = log(zinc) ~ sqrt(dist) ho una non stationary formula
# model=v.fit è i mio initial model

# Estimate the variogram from GLS residuals:
v.gls<-variogram(meuse.gstat)
v.gls.fit <- fit.variogram(v.gls, vgm(1, "Sph", 800, 1))
plot(v.gls, v.gls.fit, pch = 19)

# Update gstat object with variogram model
meuse.gstat <- gstat(id = 'zinc', formula = log(zinc) ~ sqrt(dist),
                     data = meuse, nmax = 50, model=v.gls.fit, set = list(gls=1))

## I have to define the covariate in s_0
s0.vec <- as.vector(slot(s0.new,'coords'))
# distance to the river: calculate the distance between s0 and each point of
# the river, then select the minimum
s0.dist <- min(rowSums(scale(meuse.riv,s0.vec)^2)) 
s0.new <- as.data.frame(c(s0.new,s0.dist))
names(s0.new) <- c('x','y','dist')
coordinates(s0.new) <- c('x','y')
s0.new <- as(s0.new, 'SpatialPointsDataFrame')
s0.new

# il resto come prima: il blue toglie lo shift
predict(meuse.gstat, s0.new)
predict(meuse.gstat, s0.new, BLUE = TRUE)
# prediction over the entire grid
lz.uk <- predict(meuse.gstat, meuse.grid, BLUE=FALSE)
# estimate of the mean over the entire grid
lz.uk.BLUE <- predict(meuse.gstat, meuse.grid, BLUE=TRUE)
spplot(lz.uk[,1], main = 'Universal Kriging, gstat')
spplot(lz.uk.BLUE[,1], main = 'Universal Kriging - drift , gstat')
# se cmabiano molto il drift è importante

############################# salvavita ############################
# create dummy: 0 = urban, 1 = vegetation
DUMMY <- rep(0,length(D))
DUMMY[which(D=='V')] <- 1
data <- data.frame(cbind(Bq,Long,Lat,DUMMY))
names(data) <- c('Bq','Long','Lat','D')
coordinates(data) <- c('Long','Lat')

# coefficient of the linear model: 
# it sufficies to estimate the drift at two locations where we have observations,
# with D=U and D=V
# data[1,] = urbane
# data[6,] = vegetation
g.tr <- gstat(formula = Bq ~ D, data = data, model = v.fit1)
predict(g.tr, data[1,], BLUE = TRUE)
predict(g.tr, data[6,], BLUE = TRUE)

# costruire un nuovo dato
D.s0=0.1970
s0=as.data.frame(matrix(c(0.3,0.24,D.s0),1,3))
names(s0)=c('X','Y','D')
coordinates(s0)=c('X','Y')


#### Functional Analisys: Smoothing ############################
# Upload noisy data
noisycurve <- read.table("noisycurvebis.txt",header=T)
Xobs0 <- noisycurve$X0  # observation cioe la y della fun
abscissa <- noisycurve$Abscissa  #la x della fun
NT <- length(abscissa) # number of locations of observations, quante x ho

x11()
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
plot(abscissa,Xobs0,xlab="t",ylab="observed data", type = "l")

rappincX1 <- (Xobs0[3:NT]-Xobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])
rappincX2 <- ((Xobs0[3:NT]-Xobs0[2:(NT-1)])/(abscissa[3:NT]-abscissa[2:(NT-1)])-(Xobs0[2:(NT-1)]-Xobs0[1:(NT-2)])/(abscissa[2:(NT-1)]-abscissa[1:(NT-2)]))*2/(abscissa[3:(NT)]-abscissa[1:(NT-2)])


# Set parameters
library(fda)
m <- 5           # spline order 
degree <- m-1    # spline degree 
nbasis <- 9      # ne provo uno a caso, in generale fai CV
# Create the basis
basis <- create.bspline.basis(rangeval=c(0,1), nbasis=nbasis, norder=m)

Xsp <- smooth.basis(argvals=abscissa, y=Xobs0, fdParobj=basis)
Xsp0 <- eval.fd(abscissa, Xsp$fd) #  the curve smoothing the data
Xsp1 <- eval.fd(abscissa, Xsp$fd, Lfd=1) # first derivative
Xsp2 <- eval.fd(abscissa, Xsp$fd, Lfd=2) # second derivative
df <- Xsp$df   #  the degrees of freedom in the smoothing curve  

x11(width = 14)
par(mfrow=c(1,3))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
legend("topleft", legend = c("noisy data","estimated curve"), col = c("black","blue"), lwd = c(1,3,2))
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsp2 ,type="l",col="blue",lwd=2)

############################# Approximate pointwise confidence intervals ####
# As in linear models, we can estimate the variance of x(t) as
# sigma^2*diag[phi*(phi'phi)^{-1}(phi)']
basismat <- eval.basis(abscissa, basis)  # ti dice quanto vale una base in ogni punto
S <- basismat%*%solve(t(basismat)%*%basismat)%*%t(basismat) #projection operator 
sum(diag(S)) #fa 9
sigmahat <- sqrt(sum((Xsp0-Xobs0)^2)/(NT-df)) #estimate of sigma
lb <- Xsp0-qnorm(0.975)*sigmahat*sqrt(diag(S))
ub <- Xsp0+qnorm(0.975)*sigmahat*sqrt(diag(S))

x11()
plot(abscissa,Xsp0,type="l",col="blue",lwd=2,ylab="")
points(abscissa,lb,type="l",col="blue",lty="dashed")
points(abscissa,ub,type="l",col="blue",lty="dashed")
# ricorda che sono pointwise! per qeusto sono puntini , non formano una curva

# Undersmoothing: number of basis too high, o viceversa
# generalized cross-validation, come scegliere l'ottimo nbasis
nbasis <- 6:30
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(0,1), nbasis[i], m)
  gcv[i] <- smooth.basis(abscissa, Xobs0, basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)] #ti dice quale nbasis usare
abline(v = nbasis[which.min(gcv)], col = 2)

### Bias-Variance tradeoff ###
# usalo se hai la truecurve, improbabile all'esame
# riga 243 fda smoothing lab 14

############################# SMOOTHING SPLINES #### 
#more flexibility della regression
breaks <- abscissa
basis <- create.bspline.basis(breaks, norder=m)
functionalPar <- fdPar(fdobj=basis, Lfdobj=3, lambda=1e-8)  
# functional parameter, having arguments: 
# basis, order of the derivative to be penalized, smoothing parameter.
Xss <- smooth.basis(abscissa, Xobs0, functionalPar)
Xss0 <- eval.fd(abscissa, Xss$fd, Lfd=0)
Xss1 <- eval.fd(abscissa, Xss$fd, Lfd=1)
Xss2 <- eval.fd(abscissa, Xss$fd, Lfd=2)
df <- Xss$df   #  the degrees of freedom in the smoothing curve
gcv <- Xss$gcv  #  the value of the gcv statistic

# generalized cross-validation (lambda piccolo undersmoothing, grande over)
lambda <- 10^seq(-12,-5,by = 0.5)  # prendi i punti equispaziati in scala log
gcv <- numeric(length(lambda))
for (i in 1:length(lambda)){
  functionalPar <- fdPar(fdobj=basis, Lfdobj=3, lambda=lambda[i])  
  gcv[i] <- smooth.basis(abscissa, Xobs0, functionalPar)$gcv
}
par(mfrow=c(1,1))
plot(log10(lambda),gcv)
lambda[which.min(gcv)] #lui è l'ottimo lambda, rifai tutto con lui

# i plot abbastanza intuitivi sono funzioni: fai i raw data e l'approx
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xss0 ,type="l",col="blue")

############################# LOCAL POLYNOMIAL REGRESSION #################################
library(KernSmooth)
m <- 5           # order of the polynomial
degree <- m-1    # degree of the polynomial
bw <- 0.05 # bandwidth
Xsm0 <- locpoly(abscissa, Xobs0, degree=degree,
                bandwidth=bw, gridsize=length(abscissa), 
                range.x=range(abscissa))
# locpoly fa local polinomial regression
Xsm0 <- Xsm0$y  # ha le ascissa che gli passo e la stima delle y

# Estimate of the derivatives, basta mettere drv=quello che vuoi
Xsm1 <- locpoly(abscissa,Xobs0,drv=1,degree=degree,bandwidth=bw,
                gridsize=length(abscissa), range.x=range(abscissa))
Xsm1 <- Xsm1$y

Xsm2 <- locpoly(abscissa,Xobs0,drv=2,degree=degree,bandwidth=bw,
                gridsize=length(abscissa), range.x=range(abscissa))
Xsm2 <- Xsm2$y

#plot della derivata
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsm1 ,type="l",col="blue",lwd=2)

# Changing the bandwidth: larger value , i will be more smooth
# a too small bandwidth leads to undersmoothing

############################# SMOOTHING for monotone curves ####
ageRng <- range(age)
norder <- 6
nbasis <- nage - 2 + norder 
wbasis <- create.bspline.basis(rangeval = ageRng, nbasis = nbasis, 
                               norder = norder, breaks = age)
# breaks is a vector specifying the break points defining the b-spline.
# in pratica gli dico mettimi i nodi in quei punti li, guarda l help

# We construct the functional parameter with penalization of the third 
# derivative

Lfdobj <- 3          
lambda <- 10^(-0.5)  
cvecf <- matrix(0, nbasis, ncasef) # this is used as initial value 
# for the numerical techniques
Wfd0 <- fd(coef = cvecf, basisobj = wbasis)
growfdPar <- fdPar(fdobj = Wfd0, Lfdobj = Lfdobj, lambda = lambda)
##
growthMon <- smooth.monotone(argvals = age, y = hgtf, WfdParobj = growfdPar)

Wfd <- growthMon$Wfd
betaf <- growthMon$beta
hgtfhatfd <- growthMon$yhatfd

velocfdUN <- deriv.fd(expr = hgtfhatfd, Lfdobj = 1)   # utile per fare derivata di
velocmeanfdUN <- mean.fd(velocfdUN)  # un functional object

accelfdUN <- deriv.fd(expr = hgtfhatfd, Lfdobj = 2)
accelmeanfdUN <- mean.fd(accelfdUN)

# Extension to multidimensional curves, secondo me inutile

############################# K-mean alignment ###################################
library(fdakma)
# Plot of original functions, t fa il trasposto
x11()
matplot(t(x),t(y0), type='l', xlab='x', ylab='orig.func')
title ('Original functions')

# Without alignment, let's try with 3 clusters
set.seed(4)
fdakma_example_noalign_0der <- kma(
  x=x, y0=y0, n.clust = 3, 
  warping.method = 'NOalignment', 
  similarity.method = 'd0.pearson',   # similarity computed as the cosine
  # between the original curves 
  # (correlation)
  center.method = 'k-means'
  #,seeds = c(1,11,21) # you can give a little help to the algorithm...
)

kma.show.results(fdakma_example_noalign_0der)
fdakma_example_noalign_0der$labels

# i comandi sono sempre questi, cambiano solo i parametri:
# - d1.perason da una similarity matrix che tiene conto delle derivate

table(fdakma_example_noalign_0der$labels,
      fdakma_example_noalign_1der$labels, dnn=c("0der", "1der"))
# questa table confronta 2 clusterizzazioni

# alignment
set.seed(1)
fdakma_example <- kma(
  x=x, y0=y0, y1=y1, n.clust = 2, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',  # similarity computed as the cosine
  # between the first derivatives 
  # (correlation)
  center.method = 'k-means'
  #seeds = c(1,21) # you can give a little help to the algorithm...
)

kma.show.results(fdakma_example)
fdakma_example$labels

# Total shifts and dilations applied to the original 
# abscissa to obtain the aligned abscissa
fdakma_example$shift  #distanza dal centro del cluster
fdakma_example$dilation

# How to choose the number of clusters and the warping method, ci mette 2 milioni di anni
kma.compare_example_3 <- kma.compare (
  x=x, y0=y0, y1=y1, n.clust = 1:3, 
  warping.method = c("NOalignment", "shift", "dilation", "affine"), 
  similarity.method = 'd1.pearson',
  center.method = 'k-means', 
  seeds = c(1,21,30), # numero di curva che usi come rapresentative del cluster
  plot.graph=TRUE)    # cioe curva 1 rappresenta il primo cluster, curva 21 il secondo ecc
# questa non solo sceglie il numero di cluster ma anche il metodo: vince affine
# without alignment molto peggio

############################# Functional principal component analysis (and fourier) ####
library(fda)
# il data_w deve essere fatto con: p covariate n funzioni (I GIORNI SONO LE FUNZIONI)
matplot(data_W,type='l',main='Canadian temperature',xlab='Day',ylab='Temperature')

# First of all we smooth the data. We choose a Fourier basis
time <- 1:365
# high dimensional basis (interpolating/overfitting)
# low dimensional basis (loss of information)
basis.3 <- create.fourier.basis(rangeval=c(0,365),nbasis=109)
data_W.fd.3 <- Data2fd(y = data_W,argvals = time,basisobj = basis.3)
plot.fd(data_W.fd.3)

# i primi 3 coef delle prime 2 covariate
names(data_W.fd.3)
data_W.fd.3$coefs[1:3,1:2]

library(fields)
plot.fd(data_W.fd.3)  # mean
lines(mean.fd(data_W.fd.3),lwd=2)
eval.3 <- eval.fd(time,data_W.fd.3)   # covariance
image.plot(time,time,(cov(t(eval.3))[1:365,]))
# lo smoothing influenza la pca, cioe la covariance

### new dataset
time <- seq(0,350,by=7) # se aggiungi lui nell asse x non hai piu
# il numero della covariata, ma il tempo
matplot(data_L,type='l',main='Lip data',ylab='Position',
        xlab='Time (millisec.)')

# Since the Fourier basis creates periodical data, it modifies the
# data near the boundaries of the domain, to make them periodical
# Better to use a b-spline basis (cambia solo il nome del comando)
basis <- create.bspline.basis(rangeval=c(0,350),nbasis=21)
data_L.fd <- Data2fd(y = data_L,argvals = time,basisobj = basis)
plot.fd(data_L.fd, main="B-splines")

### pca vera sul primo dataset
# interpolated data (Choice 1), in alternativa usa data_W.fd.2, smooth data
pca_W.1 <- pca.fd(data_W.fd.1,nharm=5,centerfns=TRUE) # con true centri i dati
# proporzione di varianza spiegata dalle prime 5 princ comp
names(pca_W.1)
pca_W.1$varprop
# valore della varianza delle prime 5 princ comp
pca_W.1$values[1:5]
# scree plot
# pca.fd computes all the 365 eigenvalues, but only the first N-1=34 are non-null
plot(pca_W.1$values[1:35],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_W.1$values)[1:35]/sum(pca_W.1$values),xlab='j',ylab='CPV',ylim=c(0.8,1))

# first two FPCs, equivalente al barplot dei loadings, ma ora i loading sono functions!
# interpretabili co,me i vecchi loadings se hanno una shape particolare
par(mfrow=c(1,2))
plot(pca_W.1, nx=100, pointplot=TRUE, harm=c(1,2), expand=0, cycle=FALSE)

# scatter plot of the scores (SE C'è UNA CAT AGGIUNGI COL=CAT)
par(mfrow=c(1,2))
plot(pca_W.1$scores[,1],pca_W.1$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)

plot(pca_W.1$scores[,1],pca_W.1$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2",xlim=c(-400,250))
text(pca_W.1$scores[,1],pca_W.1$scores[,2],dimnames(data_W)[[2]], cex=1)

# outlier: Resolute (35)
layout(1)
matplot(eval.1,type='l')
lines(eval.1[,35],lwd=4, col=2) #temperature profile for Resolute, outlier

coord <- CanadianWeather$coordinates
coord[,2] <- -coord[,2]
plot(coord[,2:1],col=0)
text(coord[,2:1],rownames(coord))
