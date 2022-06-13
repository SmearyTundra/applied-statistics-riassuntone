horsecolic <- read.table("horsecolic.txt", header = T)
horsecolic_pain <- horsecolic[which(horsecolic$Pain=="Yes"),1:4]
horsecolic_no_pain <- horsecolic[which(horsecolic$Pain=="No"),1:4]

np <- dim(horsecolic_pain)[1]
nn <- dim(horsecolic_no_pain)[1]
p <- dim(horsecolic_pain)[2]
n <- np + nn


load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
library(mvnormtest)
mcshapiro.test(horsecolic_no_pain)
mcshapiro.test(horsecolic_pain)
#we can assume gaussianity amongst the two groups
S_pain <- cov(horsecolic_pain)
S_no_pain <- cov(horsecolic_no_pain)
S_pooled <- ((np - 1)*S_pain + (nn - 1)*S_no_pain)/(n - 2)

x11()
image(S_pain)
x11()
image(S_no_pain)
x11()
image(S_pooled)  

sum(diag(S_pain))
sum(diag(S_no_pain))
sum(diag(S_pooled))
#homoschedasticity is plausible

alpha_bonf <- 0.01/4

#we have an exact formula for those intervals:
#a'(smean1 - smean2) ~ N(a'(mean1 - mean2), a'(1/n1 + 1/n2)SIGMA a)
#-> a'(smean1 - smean2)/sqrt(a'(1/n1 + 1/n2)SIGMA a) ~ N(0, I)
#Pivotal stat:
#a'(smean1 - smean2 - (deltamean))/sqrt(a'(1/n1 + 1/n2)Spooled a) ~ t(n - g)  (Spooled has n - 2 degrees of freedom!, not n-1!!!)

smean_pain <- sapply(horsecolic_pain, mean)
smean_no_pain <- sapply(horsecolic_no_pain, mean)

cfrt <- qt(1 - alpha_bonf/2, n - 2)

IC <- cbind(inf = smean_pain - smean_no_pain - cfrt*sqrt(diag(S_pooled)*(1/np + 1/nn)),center = smean_pain - smean_no_pain ,sup = smean_pain - smean_no_pain + cfrt*sqrt(diag(S_pooled)*(1/np + 1/nn)))
IC
#So there is a significant difference in terms of pulse and respiratory rate

library(MASS)

attach(horsecolic)
#we can use lda since we've already verified gaussianity and the same covariance structure amongst the two groups
#we are assuming prior probabilities to be the empirical ones and the costs of misclassification to be the same
ld <- lda(Pain ~ Pulse + Respiratory.rate)
ld
#from this we get the estimated priors: NO = 73%
#also the groups means, obtainable also by smean_pain and smean_no_pain

x  <- seq(min(Pulse), max(Pulse), length=200)
y  <- seq(min(Respiratory.rate), max(Respiratory.rate), length=200)
xy <- expand.grid(Pulse=x, Respiratory.rate=y)

z  <- predict(ld, xy)$post

x11()
plot(Pulse, Respiratory.rate)
points(horsecolic_no_pain$Pulse,horsecolic_no_pain$Respiratory.rate, col = "blue")
points(horsecolic_pain$Pulse,horsecolic_pain$Respiratory.rate, col = "red")
z1 <- z[,1] - z[,2]  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- z[,2] - z[,1]  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  


ldpred <- predict(ld, horsecolic)
errors <- (ldpred$class != Pain)
sum(errors)
APER   <- sum(errors)/n
APER


#we could use qda, with the same assumptions
qd <- qda(Pain ~ Pulse + Respiratory.rate)
qd

zbis  <- predict(qd, xy)$post

x11()
plot(Pulse, Respiratory.rate)
points(horsecolic_no_pain$Pulse,horsecolic_no_pain$Respiratory.rate, col = "blue")
points(horsecolic_pain$Pulse,horsecolic_pain$Respiratory.rate, col = "red")
zbis1 <- zbis[,1] - zbis[,2]  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- zbis[,2] - zbis[,1]  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    
contour(x, y, matrix(zbis1, 200), levels=0, drawlabels=F, add=T)  


qdpred <- predict(qd, horsecolic)
errors <- (qdpred$class != Pain)
sum(errors)
APERq   <- sum(errors)/n
APERq



#analogous classifications, in order to select one should study their robustness