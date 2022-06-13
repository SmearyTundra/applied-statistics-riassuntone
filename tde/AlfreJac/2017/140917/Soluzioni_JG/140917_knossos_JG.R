knossos <- read.table("knossos.txt", header = T)
attach(knossos)
summary(X)
summary(Y)

#compute dissimilarity matrix
diss <- dist(knossos, method = 'euclidean')
cl <- hclust(diss, method = 'complete')
x11()
plot(cl, main='Euclidian - complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')

clust <- cutree(cl, k = 2)
mean_1 <- sapply(knossos[which(clust == 1), ], mean)
mean_2 <- sapply(knossos[which(clust == 2), ], mean)
n_1 <- length(which(clust == 1))
n_2 <- length(which(clust == 2))
cop <- cophenetic(cl)
cop_coeff <- cor(diss, cop)
#0.87 -> pretty good

#The test we want to perform is H0: mean position(i.e. site position) for the two clusters is the same vs H1
#It's a MANOVA
#we first of all need gaussianity
an <- manova(as.matrix(knossos) ~ clust)
summary.aov(an)
summary(an)
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
library(mvnormtest)
mcshapiro.test(an$residuals)
mcshapiro.test(an$residuals[which(clust==1),])
mcshapiro.test(an$residuals[which(clust==2),])
cov1 <- cov(an$residuals[which(clust==1),])
cov2 <- cov(an$residuals[which(clust==2),])
#we can assume gaussianity in any way, and homoschedasticity holds
               