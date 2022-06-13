# EXAM 14/9/2017 ex 3

setwd("D:/Dati/Ele/poli/anno4/applied statistics/OldExams/tde")
load("D:/Dati/Ele/poli/anno4/applied statistics/R lab/LAB_5/mcshapiro.test.RData")

library(mvtnorm)
library(rgl)
library(carData)
library(car)

# Recent studies in the Knossos area have discovered numerous fragments of
# amphoras. The positions of these fragments (file knossos.txt) with respect to
# the centre of Knossos Palace (considered as origin of the coordinate system),
# suggest the existence of a second site of archeological interest.

knossos <- read.table('knossos.txt', header=T)
head(knossos)

# a) Identify two clusters of locations through a hierarchical clustering method
#    (Euclidean distance and complete linkage). Report the estimate of the mean
#    within the groups, their size, and compute the cophenetic coefficient.

# Hierarchical clustering to identify the 2 clusters
data.e <- dist(knossos, method = 'euclidean')
data.ec <- hclust(data.e, method = 'complete')
cluster.ec <- cutree(data.ec, k=2)

x11()
par(mfrow=c(1,2))
plot(data.ec, main='Euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(data.ec, k=2)
plot(knossos, col=as.vector(cluster.ec)+1, pch=19)

knossos1 <- knossos[which(cluster.ec==1),]
knossos2 <- knossos[which(cluster.ec==2),]

# Mean within the 2 groups
m1 <- sapply(knossos1, mean)
m1 # (X, Y) = (2.962927, 1.503984) 
m2 <- sapply(knossos2, mean)
m2 # (X, Y) = (0.03269231, -0.02320513)

# Size of the 2 groups
n1 <- dim(knossos1)[1]
n1 # 123
n2 <- dim(knossos2)[1]
n2 # 78

# Cophenetic coefficient
coph <- cophenetic(data.ec)
ec <- cor(data.e, coph)
ec # 0.8709168 -> Good clustering

# _________________________________________________________________________

# b) Assume the identified groups to be independent. Having introduced and
#    verified the needed assumptions, test the hypothesis according to which
#    only one archeological site exists in the Knossos area. Write a report of
#    max lines to the archeologists summarising the results of the analysis.

# TEST FOR THE MEANS OF 2 INDEPENDENT GAUSSIAN POPULATIONS

# Check the assumption of Gaussianity in each group:
mcshapiro.test(knossos1) # pvalue=0.62 -> OK
mcshapiro.test(knossos2) # pvalue=0.592 -> OK

# Test: H0: There exists only one cluster   vs   H1: !H0
# i.e.  H0: m1==m2   vs   H1: m1!=m2

S1 <- cov(knossos1)
S2 <- cov(knossos2)
Sp <- ((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)
Sp.inv <- solve(Sp)
# Check the assumption on same covariance structures
list(S1=S1, S2=S2, Spooled=Sp) # Same covariance structures

delta.0 <- c(0,0)
p <- dim(knossos1)[2]

T2 <- n1*n2/(n1+n2)*(m1-m2-delta.0)%*%Sp.inv%*%(m1-m2-delta.0)
P <- 1-pf(T2/((n1+n2-2)*p/(n1+n2-1-p)), p, n1+n2-1-p)
P
# Pvalue=0 -> Reject H0 for any alpha -> The means are significantly different
# -> There exists two different archeological sites