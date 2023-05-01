##### CLASSIFIER 

#for misclassification costs look at lab 8 line 470

library(MASS)
library(class)
library(rgl)
library(mvtnorm)
library(e1071)
library(ISLR)

x <- data.frame(Infg = 0)
predict(cyto.lda, x)$class
predict(cyto.lda, x)$posterior





###############
###### LDA (same variance)
#################
rm(list=ls())
load('mcshapiro.test.RData')
tilda = '~'

sansiro <- read.table('sansiro.txt', header = T)
head(sansiro)
sansiro$opera <- factor(sansiro$opera)
levels(sansiro$opera)
prior <- c(0.03, 0.89, 0.08)

can <- which(sansiro$opera == 'cancel')
reg <- which(sansiro$opera == 'regular')
sus <- which(sansiro$opera == 'suspend')

# verify assumptions 1) e 2): 
# 1) normality (univariate) within the groups
mcshapiro.test(sansiro[can,1:2])
mcshapiro.test(sansiro[reg,1:2])
mcshapiro.test(sansiro[sus,1:2])

# 2) equal variance (univariate)
boxM(sansiro[,1:2], sansiro$opera)
# or bartlett.test(variable-s, group)
#same covariances -> we can use linear discriminant analysis

# for the multivariate test of covariance matrixes use lab 7 qualitatively

plot(sansiro[,1:2], col = sansiro$opera)

##LDA
# Linear Discriminant Analysis (LDA)

# X|Li ~ N_p(m_i , cov)

lda.san <- lda(sansiro[,1:2], sansiro$opera, prior=prior)
lda.san
Lda.san <- predict(lda.san, sansiro[,1:2])
names(Lda.san)

#APER if no prior
Lda.san$class   # assigned classes
sansiro$opera     # true labels
table(class.true=sansiro$opera, class.assigned=Lda.san$class)

errors <- (Lda.san$class != sansiro$opera)
errors
sum(errors)
length(sansiro$opera)

APER   <- sum(errors)/length(sansiro$opera)
APER

#APER if prior
G <- length(prior)
misc <- table(class.true=species.name, class.assigned=Lda.iris$class)
APER <- 0
for(g in 1:G)
    APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior[g]  


# AER leave-one-out if no prior
errors_CV <- 0
for(i in 1:150){
    LdaCV.i <- lda(iris2[-i,], species.name[-i], prior=c(50,50,50)/150)
    errors_CV <- errors_CV + 
        as.numeric(predict(LdaCV.i,iris2[i,])$class != species.name[i])
}
errors_CV

AERCV <- errors_CV / length(species.name)
AERCV




#CLASSIFICATION REGION
plot(sansiro[,1:2], main='Sansiro', xlab='temp', ylab='prec', pch=20, col = sansiro$opera)
points(lda.san$means, pch=4,col=c('black','red','green') , lwd=2, cex=1.5)
x  <- seq(min(sansiro[,1]), max(sansiro[,1]), length=200)
y  <- seq(min(sansiro[,2]), max(sansiro[,2]), length=200)
xy <- expand.grid(temperature=x, precipitation=y)

z  <- predict(lda.san, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - pmax(z[,2], z[,3])  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- z[,2] - pmax(z[,1], z[,3])  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    
z3 <- z[,3] - pmax(z[,1], z[,2])  # P_3*f_3(x,y)-max{P_j*f_j(x,y)}

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)
contour(x, y, matrix(z3, 200), levels=0, drawlabels=F, add=T)

#NEWDATA
new <- data.frame(temperature = 28, precipitation = 12)
answer <- predict(lda.san, new)







##################
###### QDA (variances can be different)
#####################
rm(list=ls())
load('mcshapiro.test.RData')
tilda = '~'
graphics.off()

trading <- read.table('trading.txt')
head(trading)

GAIN <- which(trading$gain == 1)
NOGAIN <- which(trading$gain == 0)

# verify assumptions 1) e 2): 
# 1) normality (univariate) within the groups
mcshapiro.test(trading[GAIN, 1:2])
mcshapiro.test(trading[NOGAIN, 1:2])

# 2) equal variance (univariate)
var.test(trading[GAIN, 1:2],trading[NOGAIN, 1:2])  #should be not equal

# X|Li ~ N_p(m_i , cov_i)

qda.trad <- qda(trading[,1:2], factor(trading$gain), prior=prior)
qda.trad

##Region
plot(trading[,1:2], main='Trading', xlab='Google', ylab='Apple', pch=20, col = factor(trading$gain))
points(qda.trad$means, col=c('black','red'), pch=4, lwd=2, cex=1.5)

x  <- seq(min(trading[,1]), max(trading[,1]), length=200)
y  <- seq(min(trading[,2]), max(trading[,2]), length=200)
xy <- expand.grid(google=x, apple=y)

z  <- predict(qda.trad, xy)$post  
z1 <- z[,1] - pmax(z[,2])    
z2 <- z[,2] - pmax(z[,1])

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

##Aper  (only if no prior)
Qda.pred <- predict(qda.trad, trading[,1:2])
Qda.pred$class
trading$gain
table(class.true=trading$gain, class.assigned=Qda.pred$class)

errorsq <- (Qda.pred$class != trading$gain)
errorsq

APERq   <- sum(errorsq)/length(trading$gain)

##APER (if prior)
Qda.m <- predict(qda.m)
table(class.true=vf, class.assigned=Qda.m$class)

#            class.assigned
#class.true   true false
#      true    98     2
#     false    80    20

APER  <- (2*pt+80*pf)/100   #pt and pf are the prior
APER

# Trivial classifier; classifies always as the most likely class a priori
APER.banale <- prior[1]
APER.banale

###AperCV (only if no prior)
QdaCV.iris <- qda(iris2, species.name, CV=T)
QdaCV.iris$class
species.name
table(class.true=species.name, class.assignedCV=QdaCV.iris$class)

errorsqCV <- (QdaCV.iris$class != species.name)
errorsqCV

AERqCV   <- sum(errorsqCV)/length(species.name)
AERqCV
APERq


# L1O Cross Validation AER WITH PRIOR
QdaCV <- qda(d[,1:2], d$norm_abnorm, CV=TRUE, prior = c(abn,nor))  # specify the argument CV
table(class.true=d$norm_abnorm, class.assignedCV=QdaCV$class)

AERCV   <- (32*abn)/70 + (8*nor)/80
cat("AERCV: ", AERCV)


# ORRRR
qn <- nrow(d)
errors_CV_AB_NO <- 0
errors_CV_NO_AB <- 0
correct_NO <- 0
correct_AB <- 0
for(i in 1:n){
    qdacv <- qda(d[-i,1:2], d[-i,3], prior=c(abn,nor))
    errors_CV_AB_NO <- errors_CV_AB_NO + 
        (as.numeric(predict(qdacv,d[i,1:2])$class == "NO" && d$norm_abnorm[i] == "AB"))
    errors_CV_NO_AB <- errors_CV_NO_AB + 
        (as.numeric(predict(qdacv,d[i,1:2])$class == "AB" && d$norm_abnorm[i] == "NO"))
    correct_NO <- correct_NO + 
        (as.numeric(predict(qdacv,d[i,1:2])$class == "NO" && d$norm_abnorm[i] == "NO"))
    correct_AB <- correct_AB + 
        (as.numeric(predict(qdacv,d[i,1:2])$class == "AB" && d$norm_abnorm[i] == "AB"))
    
}

tab <- matrix(c(correct_AB, errors_CV_AB_NO, errors_CV_NO_AB, correct_NO), ncol=2, byrow=TRUE)

colnames(tab) <- c('AB','NO')
rownames(tab) <- c('AB','NO')
tab <- as.table(tab)
tab

and then look at the table





###### KNN + CV ####
library(class)
MISS <- NULL
set.seed(321)
for (k in 10:30) {
    debris.knn <- knn.cv(train = debris[,1:2], k = k, cl = debris$risk)
    misclassif.error <- length(which(debris$risk != debris.knn))
    MISS <- cbind(MISS, misclassif.error)
}

error.rate <- min(MISS)/length(db$risk)

kopt <- which.min(MISS)+9

debris.knn.plot <- knn(train = debris[,1:2], test=xy, k = kopt, cl = debris$risk)
#plot(xy, col=debris.knn.plot)
plot(debris[,1:2], col = ifelse(debris$risk=='L', 'red', 'black'), pch = 19)
z  <- as.numeric(debris.knn.plot)
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5), drawlabels=F, add=T)

new.datum <- data.frame(x=1, y=-4)
pred <- knn(train = debris[,1:2], test=new.datum, k = kopt, cl = debris$risk)






SVM

# Fit a Support Vector Machine (kernel = "radial") given a cost C
svmfit <- svm(y~., data=dat [train ,], kernel ='radial', gamma =1, cost =1)
summary(svmfit)

# Plot the SVM
plot(svmfit , dat, col =c('salmon', 'light blue'), pch=19, asp=1)