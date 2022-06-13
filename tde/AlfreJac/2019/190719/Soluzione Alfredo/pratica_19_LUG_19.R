### exam
##
setwd("Documents/Polimi/APPSTATS/")

pines <- read.table("pines.txt")
head(pines)
# test plausible mean
"H0: mu = c(14.2350, 42.4520). H1: mu != that value"
"Reject H0 if T2 higher than the soglia, knowing it is an F distribution, with
99% probability of stating they are different when they actually are not"
Xb <- sapply(pines, mean)
S <- cov(pines)
n <- dim(pines)[1]
p <- dim(pines)[2]
n
p
mu <- c(14.2350, 42.4520)

T2 <- n*t(Xb - mu) %*% solve(S) %*% (Xb - mu)
# obtain the soglia: confer the F distribution 99th percentile
cfr.F <- (n-1)*p/(n-p) * qf(.99, p, n-p)
T2 < cfr.F
"We assume the proposed mu as the mean of the population, given that with the
data we have no statistical evidence to reject it."
"In fact, the p value, i.e. the probability of observing such value if the H0
is true, is reasonable according to the set standard"
pval <- 1-pf(T2, p, n-p, lower.tail=T)
pval
"NOw verify the assumptions: w
  every sample (row) of the dataset is an iid sample CHECK
  it is a multivariate normal population: we test it with the monte carlo
shapiro test: even the 'worst' direction in terms of W statistic a normal distribution
since all lineaer combinations of a normal distribution have to be normal"
mcshapiro.test(pines)
"at test level of 0.01 we do not reject it, so we say it is multivariate normal"

#### b####
" this is a prediction interval:the region in which 99% of the points of the
normal distribution will be"
"From a result we know the mahalanobis distance from a sample to its mean is 
chi quadro(p), so
We also assume the true mean and the true covariance matrix the ones estimated
by LLN"
raggio <- qchisq(0.99, p, lower.tail = T)
# direction of axes
eigen(S)$vectors[,1]
eigen(S)$vectors[,2]
# length of axes
eigen(S)$values[1] * sqrt(raggio)
eigen(S)$values[2] * sqrt(raggio)

# qualitative plot of the region
library(car)
x11()
plot(pines, asp = 1)
ellipse(mu, shape = S,
        sqrt(raggio), col = 'blue', lty = 2, center.pch = 16)

###############
### Problem 2
#######
buoys <- read.table("buoys.txt")
# compute the ditances
buoys.e <- dist(buoys[,c(1,2)], method='euclidean')
# now perform hierarchical clustering using ward's method (minimizing sum of squares)
# within group at each agglomeration
buoys.ew <- hclust(buoys.e, method='ward.D2') 
x11()
plot(buoys.ew, hang=-0.1, labels=FALSE, main='ward', xlab='', sub='')
rect.hclust(buoys.ew, k=2)
rect.hclust(buoys.ew, k=3)
"I deem appropiate 2 the number of clusters because it is the k at which the
change in distance takes a lot to change, i.e. the gomito on the n clusters
vs distance plot"
## DETERMING OPTIMAL CLUSTERS
# Elbow Method
fviz_nbclust(buoys[,c(1,2)], FUN = hcut, method = "wss")


clusters <- cutree(tree = buoys.ew, k=2)
# numerosity
table(clusters)
n1 <- table(clusters)[1]
n2 <- table(clusters)[2]
which(clusters==1)
mean1 <- sapply(buoys[which(clusters==1),c(1,2)], mean)
mean1
mean2 <- sapply(buoys[-which(clusters==1),c(1,2)], mean)
mean2

## b)

shapiro.test(buoys[which(clusters==1),3])
shapiro.test(buoys[which(clusters==2),3])
bartlett.test(buoys$DO,  clusters)
t.test(buoys$DO, clusters, var.equal = T)

mcshapiro.test(fit$residuals)
"We perform a MANOVA, we need equal covariances, multivariate populations,
iid samples"
"H0:"
"In the case of 1 factor with two levels, we can do it 'by hand'"
# verify gaussianity on me
mcshapiro.test(buoys[which(clusters==1),c(1,2)])
mcshapiro.test(buoys[-which(clusters==1),c(1,2)])
  
############
# BONUS
##########
# avevo capito male e quindi ho fatto un'altra cosa
# noe equal covariances
library(biotools)
buoys$cl <- 1:(n1+n2) %in% which(clusters==1)
boxM(buoys[,c(1,2)], as.factor(buoys$cl))
"However, we know that even if we reject box's M test, the MANOVA is not
damaged, especially with equal sample sizes.
Also, usually the rule of thumb is used: of the discrempacy of a multiple of 
4 of the other, it is too much"
S1 <- cov(buoys[which(clusters==1),c(1,2)])
S2 <- cov(buoys[-which(clusters==1),c(1,2)])
S1 * 1000
S2 * 1000
# we cannot say the are the same covariance matrix, so we pass
# the argument that the covariance matrices are different, 
# so we use the solution to the behrens fisher problem proposed in the book
# which is used by manva
fit <- manova(as.matrix(buoys[,c(1,2)]) ~  buoys$cl)
summary.manova(fit,test="Wilks")
mcshapiro.test(fit$residuals)





#############
### PROBLEM 3
#############
piad <- read.table("piadeina.txt")
head(piad)
names(piad)
# make dummy
piad$Day.of.Week <- as.factor(piad$Day.of.Week )

attach(piad)
fm <- lm(Sales ~., data = piad)
summary(fm) 

fm2 <- lm(Sales ~ Day.of.Week + Bread.Sold + Sandwich.Sold + Focaccia.Sold 
         + Piadina.Sold + Juices.Sold +
           Total.Soda.and.Coffee.Sold + Max.Daily.Temperature, data = piad)
summary(fm2)
fm3 <- lm(Sales ~ Day.of.Week + Bread.Sold + Sandwich.Sold + Focaccia.Sold 
          + Piadina.Sold + Juices.Sold +
            Total.Soda.and.Coffee.Sold, data = piad)
summary(fm3)
fm3 <- lm(Sales ~ Day.of.Week + Sandwich.Sold + Focaccia.Sold 
          + Piadina.Sold + Juices.Sold +
            Total.Soda.and.Coffee.Sold, data = piad)
summary(fm3)
# we see we will need a whole algorithm so we stop here.
# we have to verify homoschedasticity and 
# standard deviation
sqrt((t(fm3$residuals) %*% fm3$residuals)/fm3$df.residual)
x11()
par(mfrow=c(2,2))
plot(fm)
"So we see approximate homoschedasticity, no pattern of residuals nor
studentized residuals vs order, no significant leverages st that cook's distance
is too much"
shapiro.test(fm3$residuals)
" Residuals normal at statistical level alpha = 0.01"

# b) lasso
piad$Fri <- as.integer(piad$Day.of.Week == "Fri")
piad$Mon <- as.integer(piad$Day.of.Week == "Mon")
piad$Tue <- as.integer(piad$Day.of.Week == "Tue")
piad$Wed <- as.integer(piad$Day.of.Week == "Wed")
library(glmnet)
attach(piad)
model.matrix(Sales ~ ., data = piad)
piad <- piad[,-1]
lasso.mod <- glmnet(x = model.matrix(Sales ~ ., data = piad), y = piad$Sales,alpha=1,lambda=5)
lasso.mod$beta

# c) Cross validation
"ten fold cross validation"
lambda.grid <- seq(0, 100, 2)
cv.lasso <- cv.glmnet(model.matrix(Sales ~ ., data = piad),
                       piad$Sales, alpha = 1, lambda=lambda.grid,nfolds = 5)
bestlam.lasso <- cv.lasso$lambda.1se
bestlam.lasso
x11()
plot(cv.lasso)
coef.lasso <- predict(lasso.mod, s=bestlam.lasso, type = 'coefficients')
coef.lasso 

#############
# Problem 4
################

water <- read.table("watertemp.txt")
head(water)
library(fda)
nbasis <- 45
fbasis <- create.fourier.basis(nbasis = 45, rangeval = c(0,365))
st.basis <- eval.basis(evalarg = 1:365, fbasis)

x11()
plot(st.basis)
# transpose water as every iid unit here is a "column"
water.fd <- Data2fd(y = t(water[,-366]), argvals = 1:365,basisobj = fbasis)
x11()
plot(water.fd)
water.fd$coefs[1:3, 1:2]


###
## b) 
###
library(fields)
pca.of.smoothed <- pca.fd(water.fd, nharm = 5)
pca.of.smoothed$values
cumsum(pca.of.smoothed$values)/sum(pca.of.smoothed$values)
x11()
plot(pca.of.smoothed$values,xlab='j',ylab='Eigenvalues')
plot(cumsum(pca.of.smoothed$values)/sum(pca.of.smoothed$values),xlab='j',ylab='CPV',ylim=c(0.8,1))


par(mfrow=c(1,3))
plot.pca.fd(pca.of.smoothed, nx=100, pointplot=TRUE,
            harm=c(1,2,3), expand=0, cycle=FALSE)
# 1st component
"Stations that are on a warmer zone, at lower latitudes on northern hemisphere,
(so temperature rises), and also summer comes a little later, so could be 
something that is also implied by the latitude"
# 2nd component
"Earlier summer, could be the effect of a sea current, and then a 
a much less cold winter, could be points closer perhaps to continet"
# 3rd component
"Milder places with less hot summer and less cold winter"
### c
### take into account the zone
colors <- ifelse(water[,366]=="Deep", "Red",
                 ifelse(water[,366]=="Medium", "Blue", "Pink"))
x11()
plot(pca.of.smoothed$scores[,1],pca.of.smoothed$scors[,2],
     xlab="Scores FPC1",ylab="Scores FPC2",lwd=2, col=colors)
legend(5, 5, legend=unique(water[,366]), col = c("Red", "Blue", "Pink"))
"So at deeper waters, we see less hot summers and more cold winters, 
and also later summer, could make sense as the changes until the and "
"So pretty much a linear behavior on the scores of the PC as the stations
are less deep"

