geisha <- geisha[,-1]#delete the tag
summary(geisha$duration)
summary(geisha$time)
#variables are comparable, no need to standardize
diss <- dist(geisha)#euclidian method implicit
n <- dim(geisha)[1]
p <- dim(geisha)[2]
x11()
image(1:n, 1:n, as.matrix(diss))

link <- hclust(diss, method = 'single')
x11()
plot(link,main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
clust <- cutree(link, k = 2)
size1 <- length(which(clust == 1))
size2 <- length(which(clust == 2))
mean1 <- colMeans(geisha[which(clust == 1),])
mean2 <- colMeans(geisha[which(clust == 2),])

coph <- cophenetic(link)
x11()
plot(geisha)
points(geisha[which(clust==1),], col = 'red')
points(geisha[which(clust==2),], col = 'blue')

#it is a very unstable clustering, as it is clear from the cophenetic distances, moreover given the more regular shapes of the data clouds, average or complete linkage might better capture the clusters

link_a <- hclust(diss, method = 'average')
x11()
plot(link_a,main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')

link_c <- hclust(diss, method = 'complete')
x11()
plot(link_c,main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')


clust_a <- cutree(link_a, k = 2)
clust_c <- cutree(link_c, k = 2)

x11()
plot(geisha)
points(geisha[which(clust_a==1),], col = 'red')
points(geisha[which(clust_a==2),], col = 'blue')

x11()
plot(geisha)
points(geisha[which(clust_c==1),], col = 'red')
points(geisha[which(clust_c==2),], col = 'blue')

#average and complete linkage yield the same result, which is indeed very robust and captures the structure of the data

size1 <- length(which(clust_a == 1))
size2 <- length(which(clust_a == 2))
#cluster 1 is the one of the successes

#we need gaussianity assumptions for the two groups
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
library(mvnormtest)
mcshapiro.test(geisha[which(clust_c==1),])
mcshapiro.test(geisha[which(clust_c==2),])
#we can assume that the two groups are generated from two gaussian models
#we need to verify wether they have the same variance
S1 <- cov(geisha[which(clust_c==1),])
S2 <- cov(geisha[which(clust_c==2),])
S1
S2
var.test(geisha$duration[which(clust_c==1)],geisha$duration[which(clust_c==2)])
var.test(geisha$time[which(clust_c==1)],geisha$time[which(clust_c==2)])
bartlett.test(geisha$time ~ clust_c)
x11()
image(1:2,1:2, S1)
x11()
image(1:2,1:2, S2)

# Homoschedasticity assumption does not hold
#Well, yes... but we do not really care

mean_1 <- colMeans(geisha[which(clust_c==1),])
mean_2 <- colMeans(geisha[which(clust_c==2),])
S_pooled <- ((size1-1)*S1 + (size2-1)*S2)/(size1+size2-2)


#we provide two pairwise simultaneous confidence intervals (one for the difference in means and one for the mean of group 1)
alpha_Bonf <- 0.1/2

#Pivotal stat (1/n1 + 1/n2)^(-1)[(X1 - X2) - (mu1 - mu2)]'(S_pooled^-1)()[(X1 - X2) - (mu1 - mu2)] ~ (n1 + n2 -2)p/(n1 + n2 - 1 - p)F(p, n1 + n2 - 1 - p)
#T2 Bonferroni corrected intervals for the difference in mean:
cfr.fisher <- (p*(size1+size2-2)/(size1+size2-1-p))*qf(1-alpha_Bonf,p,size1+size2-1-p)

CI.delta.T2.duration <- c(mean_1[1]-mean_2[1]-sqrt(cfr.fisher*S_pooled[1,1]*(1/size1+1/size2)),mean_1[1]-mean_2[1] ,mean_1[1]-mean_2[1] + sqrt(cfr.fisher*S_pooled[1,1]*(1/size1+1/size2)))
CI.delta.T2.time <- c(mean_1[2]-mean_2[2]-sqrt(cfr.fisher*S_pooled[2,2]*(1/size1+1/size2)),mean_1[2]-mean_2[2],mean_1[2]-mean_2[2]+sqrt(cfr.fisher*S_pooled[2,2]*(1/size1+1/size2)))

#now a simple CI for the mean, again T2 intervals with the corrected alpha
#Pivotal stat: n(X - mu)'S^-1(X - mu) ~ (n-1)p/(n-p)F(p,n-p)
cfr.fisher.mean <- ((size1 -1)*p/(size1-p))*qf(1-alpha_Bonf,p,size1-p)
CI.mean.T2.duration <- c(mean_1[1]-sqrt(cfr.fisher.mean*S1[1,1]*(1/size1)),mean_1[1],mean_1[1] + sqrt(cfr.fisher.mean*S1[1,1]*(1/size1)))
CI.mean.T2.time <- c(mean_1[2]-sqrt(cfr.fisher.mean*S1[2,2]*(1/size1)),mean_1[2],mean_1[2] + sqrt(cfr.fisher.mean*S1[2,2]*(1/size1)))

CI <- t(matrix(c(CI.delta.T2.duration,CI.delta.T2.time,CI.mean.T2.duration,CI.mean.T2.time), nrow = 3, ncol = 4))
colnames(CI) <- c("inf", "point", "sup")
rownames(CI)<-c("delta duration", "delta time", "succ duration", "succ time")

#the successfull tours are those that start, on average, approximately at 16:45 (i.e 10/20 minutes later the unsuccessfull toures)
#moreover they last approximately 40/50 minutes longer than the unsuccessful ones