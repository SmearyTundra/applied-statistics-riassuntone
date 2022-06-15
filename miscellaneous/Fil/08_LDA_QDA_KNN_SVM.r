### TOPICS:
### Linear and Quadratic Discriminant Analysis
### Support Vector Machines


library(MASS)
library(e1071)




############# LDA ###################
###--------------------------
species.name <- factor(Species, labels=c('setosa','versicolor','virginica'))

g=3 
i1 <- which(species.name=='setosa')
i2 <- which(species.name=='versicolor')
i3 <- which(species.name=='virginica')

n1 <- length(i1)
n2 <- length(i2)
n3 <- length(i3)
n <- n1+n2+n3   

# plot the data
x11()
plot(iris2, main='Iris Sepal', xlab='Sepal.Length', ylab='Sepal.Width', pch=19)
points(iris2[i1,], col='red', pch=19)
points(iris2[i2,], col='green', pch=19)
points(iris2[i3,], col='blue', pch=19)
legend("topright", legend=levels(species.name), fill=c('red','green','blue'))

m <-  colMeans(iris2)
m1 <- colMeans(iris2[i1,])
m2 <- colMeans(iris2[i2,])
m3 <- colMeans(iris2[i3,])

S1 <- cov(iris2[i1,])
S2 <- cov(iris2[i2,])
S3 <- cov(iris2[i3,])
Sp  <- ((n1-1)*S1+(n2-1)*S2+(n3-1)*S3)/(n-g)


#### Assumptions:
###------------------
# 1) if L=i, X.i ~ N(mu.i, sigma.i^2), i=A,B
# 2) sigma.A=sigma.B
# 3) c(A|B)=c(B|A) (equal misclassification costs)

# verify assumptions 1) e 2): 
# 1) normality (univariate) within the groups
shapiro.test(cyto[A,1])
shapiro.test(cyto[B,1])

# 2) equal variance (univariate)
var.test(cyto[A,1],cyto[B,1])


#### LDA (univariate) 
###------------------
lda.iris <- lda(iris2, species.name)
lda.iris
# In particular:
# - coefficients of linear discriminants: versors of the canonical directions
#   [to be read column-wise]
# - proportion of trace: proportion of variance explained by the corresponding 
#   canonical direction


#### Prediction and APER 
###------------------
Lda.iris <- predict(lda.iris, iris2)
names(Lda.iris)

# 1) Compute the APER
table(class.true=species.name, class.assigned=Lda.iris$class)
errors <- (Lda.iris$class != species.name)
APER   <- sum(errors)/length(species.name)
APER

# 2) AER L1OCV
LdaCV.iris <- lda(iris2, species.name, CV=TRUE)  # specify the argument CV
table(class.true=species.name, class.assignedCV=LdaCV.iris$class)
errorsCV <- (LdaCV.iris$class != species.name)
AERCV   <- sum(errorsCV)/length(species.name)
AERCV


#### Plot partition induced by LDA (rescaled for 2 groups)
###------------------
x11()
plot(iris2, main='Iris Sepal', xlab='Sepal.Length', ylab='Sepal.Width', pch=20)
points(iris2[i1,], col='red', pch=20)
points(iris2[i2,], col='green', pch=20)
legend("topright", legend=levels(species.name), fill=c('red','green'), cex=.7)

points(lda.iris$means, pch=4,col=c('red','green') , lwd=2, cex=1.5)
x  <- seq(min(iris[,1]), max(iris[,1]), length=200)
y  <- seq(min(iris[,2]), max(iris[,2]), length=200)
xy <- expand.grid(Sepal.Length=x, Sepal.Width=y)

z  <- predict(lda.iris, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - pmax(z[,2])  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- z[,2] - pmax(z[,1])  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    

# Plot the contour line of level (levels=0) of z1, z2, z3: 
# P_i*f_i(x,y)-max{P_j*f_j(x,y)}=0 i.e., boundary between R.i and R.j 
# where j realizes the max.
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)


















############# QDA ###################
###--------------------------
#### Assumptions:
###------------------
# 1) if L=i, X.i ~ N(mu.i, sigma.i^2), i=A,B
# 2) c(A|B)=c(B|A) (equal misclassification costs)

qda.iris <- qda(iris2, species.name)
qda.iris
Qda.iris <- predict(qda.iris, iris2)

# compute the APER
table(class.true=species.name, class.assigned=Qda.iris$class)
errorsq <- (Qda.iris$class != species.name)
APERq   <- sum(errorsq)/length(species.name)
APERq
# Remark: correct only if we estimate the priors through the sample frequencies!

# Compute the estimate of the AER by leave-out-out cross-validation 
QdaCV.iris <- qda(iris2, species.name, CV=T)
table(class.true=species.name, class.assignedCV=QdaCV.iris$class)
errorsqCV <- (QdaCV.iris$class != species.name)
AERqCV   <- sum(errorsqCV)/length(species.name)

# Plot partition
x11()
plot(iris2, main='Plot', pch=20)
points(iris2[i1,], col='red', pch=20)
points(iris2[i2,], col='green', pch=20)
legend("topright", legend=levels(species.name), fill=c('red','green'), cex=.7)

points(qda.iris$means, pch=4,col=c('red','green') , lwd=2, cex=1.5)
x  <- seq(min(iris[,1]), max(iris[,1]), length=200)
y  <- seq(min(iris[,2]), max(iris[,2]), length=200)
xy <- expand.grid(Sepal.Length=x, Sepal.Width=y) #### CHANGE THIS ONE

z  <- predict(qda.iris, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - pmax(z[,2])  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- z[,2] - pmax(z[,1])  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    

# Plot the contour line of level (levels=0) of z1, z2, z3: 
# P_i*f_i(x,y)-max{P_j*f_j(x,y)}=0 i.e., boundary between R.i and R.j 
# where j realizes the max.
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)
















############## KNN ###################
###--------------------------
k <- 7

x11()
plot(iris2, main='Plot', xlab='x1', ylab='x2', pch=20)
points(iris2[i1,], col=2, pch=20)
points(iris2[i2,], col=3, pch=20)
points(iris2[i3,], col=4, pch=20)
legend("topright", legend=levels(species.name), fill=c(2,3,4))

x  <- seq(min(iris[,1]), max(iris[,1]), length=200)
y  <- seq(min(iris[,2]), max(iris[,2]), length=200)
xy <- expand.grid(Sepal.Length=x, Sepal.Width=y) ### CHANGE THIS

iris.knn <- knn(train = iris2, test = xy, cl = iris$Species, k = k)
z  <- as.numeric(iris.knn)
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5), drawlabels=F, add=T)





















################# WITH PRIORS/MISCLASSIFICATION COSTS ###################
###------------------------------

true <- read.table('moneytrue.txt',header=TRUE)
false <- read.table('moneyfalse.txt',header=TRUE)
banknotes <- rbind(true,false)
vf <- factor(rep(c('true','false'),each=100), levels=c('true','false'))
x11()
plot(banknotes[,1:2], main='Banknotes', xlab='V1', ylab='V2', pch=20)
points(false, col='red', pch=20)
points(true, col='blue', pch=20)
legend('bottomleft', legend=levels(vf), fill=c('blue','red'), cex=.7)


# misclassification costs
c.vf <- 10       # MISCLASSIFYING FALSE AS TRUE
c.fv <- 0.05     # MISCLASSIFYING TRUE AS FALSE

#prior probabilities
pf <- 0.001
pt <- 1-0.001

# Prior modified to account for the misclassification costs
prior.c <- c(pt*c.fv/(c.vf*pf+c.fv*pt), pf*c.vf/(c.vf*pf+c.fv*pt))
prior.c
### PRIOR(TRUE) * FRACTION(COST OF MISCLASSIFYING TRUE (AS FALSE))
### PESO DI PIù UN PRIOR (COME SE CE NE FOSSERO DI PIù) SE IL COSTO DI 
### MISCLASSIFICARE QUEL PRIOR è ALTO

qda.m <- qda(banknotes, vf, prior=prior.c)
qda.m


# APER
Qda.m <- predict(qda.m)
table(class.true=vf, class.assigned=Qda.m$class)
APER  <- (2*pt+80*pf)/100     #### (TRUE CLASSIFIED AS FALSE * PRIOR(TRUE))/NUMBER(TRUE) +
APER                          #### (FALSE CLASSIFIED AS TRUE * PRIOR(FALSE))/NUMBER(FALSE) 



















################# SVM ###################
###------------------------------
# Fit the Support Vector Classifier (kernel = "linear")
# given a cost C
dat <- data.frame(x=x, y=as.factor (y))
svmfit <- svm(y~., data=dat , kernel ='linear', cost =10, scale =FALSE )
summary(svmfit)

x11()
par(mfrow=c(1,2))
plot(svmfit , dat, col =c('salmon', 'light blue'), pch=19, asp=1)



#### Tuning cost for best mod
###------------------------------
set.seed (1)
tune.out <- tune(svm,y~.,data=dat ,kernel = 'linear',
                 ranges =list(cost=c(0.001 , 0.01, 0.1, 1,5,10,100) ))
summary(tune.out)
# Extract the best model from the result of tune
bestmod <- tune.out$best.model
summary(bestmod)
plot(bestmod , dat, col =c('salmon', 'light blue'), pch=19, asp=1)



#### Prediction for a new observation
###------------------------------
table(true=data[,"y"], pred=predict(svmfit,newdata=data))

testdat <- data.frame(x1=,x2=)
ypred <- predict(bestmod,testdat)
