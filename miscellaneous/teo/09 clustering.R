#######
# CLUSTERING
######

n <- dim(data)[1] # number of samples
misc <- sample(n) # random permutation

# d = distance
# l = linkage

# distance matrix
dist.matrix <- dist(data, method='euclidean-canberra-manhattan')

# plot dist matrix
image(1:n, 1:n, as.matrix(dist.matrix), main = "metrics: NAME OF METRIC", asp = 1, xlab = "i", ylab = "j")

# clustering
cluster.dl <- hclust(dist.matrix, method = "single-average-complete-ward.D2")
names(cluster.dl) # $order

# plot dendrogram
plot(cluster.dl, main = "distance-linkage", hang = -0.1, xlab = "", labels = F, cex = 0.6, sub = "")
# with k cluster
rect.hclust(cluster.dl, k = 2-3-4)

# cut the cluster
cluster.dl.cut <- cutree(cluster.dl, k = 2)

# how many
table(cluster.dl.cut)

# plot distinction
plot(data, col = ifelse(cluster.dl.cut == 1, "red", "blue"), pch = 19) # or col = cluster.dl (+1 for changing colors)

# cophenetic matrix
coph.matrix <- cophenetic(cluster.dl)
# cophenetic coefficient
coph.coeff <- cor(dist.matrix, coph.matrix)




##### WARD + Intervals for mean and variance ####
rm(list=ls())
load('mcshapiro.test.RData')
tilda = '~'
graphics.off()

sequoia <- read.table('sequoia.txt', header=T)
head(sequoia)
plot(sequoia)

sequoia.e <- dist(sequoia, method='euclidean')
sequoia.ew <- hclust(sequoia.e, method='ward')
plot(sequoia.ew, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(sequoia.ew, k=5)
cluster.ew <- cutree(sequoia.ew, k=5)

plot(sequoia, col = cluster.ew)

#Discuss about the algorithm's goodness.
coph <- cophenetic(sequoia.ew)
cor(coph, sequoia.e) # cophenetic correlation coefficient -> the higher, the better

##
c1 <- which(cluster.ew == 1)
c2 <- which(cluster.ew == 2)
c3 <- which(cluster.ew == 3)
c4 <- which(cluster.ew == 4)
c5 <- which(cluster.ew == 5)

n1 <- length(c1) 
n2 <- length(c2) 
n3 <- length(c3) 
n4 <- length(c4) 
n5 <- length(c5) 

m1 <- colMeans(sequoia[c1,])
m2 <- colMeans(sequoia[c2,])
m3 <- colMeans(sequoia[c3,])
m4 <- colMeans(sequoia[c4,])
m5 <- colMeans(sequoia[c5,])

#bonferroni for mean and variance
shapiro.test(sequoia[c1,2])$p.value
shapiro.test(sequoia[c2,2])$p.value
shapiro.test(sequoia[c3,2])$p.value
shapiro.test(sequoia[c4,2])$p.value
shapiro.test(sequoia[c5,2])$p.value

alpha <- 0.1
k <- 10

var1 <- var(sequoia[c1,2])
cfr.t1 <- qt(1-alpha/(2*k), n1-1)
ICm1 <- c(m1[2] - sqrt(var1/n1)*cfr.t1, m1[2], m1[2]  + sqrt(var1/n1)*cfr.t1)
ICv1 <- c((n1-1)*var1/qchisq(1-alpha/k, n1-1), var1, (n1-1)*var1/qchisq(alpha/k, n1-1))

var2 <- var(sequoia[c2,2])
cfr.t2 <- qt(1-alpha/(2*k), n2-1)
ICm2 <- c(m2[2] - sqrt(var2/n2)*cfr.t2, m2[2], m2[2]  + sqrt(var2/n2)*cfr.t2)
ICv2 <- c((n2-1)*var2/qchisq(1-alpha/k, n2-1), var2, (n2-1)*var2/qchisq(alpha/k, n2-1))

# FASTER


alpha <- 0.05
k <- 6 # number of CI
g <- 3 # number of clusters

IC={}
Ps={}
for(i in 1:g){
    X <- data[cluster.dl.cut==i,1] # i need only the major axis
    n <- length(X)
    Ps <- c(Ps,shapiro.test(X)$p)
    x.mean   <- mean(X)
    x.cov    <- var(X)
    
    ICmean <- c(inf    = x.mean - sqrt(x.cov/n) * qt(1 - alpha/(2*k), n-1),
                center = x.mean,
                sup    = x.mean + sqrt(x.cov/n) * qt(1 - alpha/(2*k), n-1))
    
    ICvar <- c(inf     = x.cov*(n-1) / qchisq(1 - alpha/(2*k), n-1),
               center  = x.cov,
               sup     = x.cov*(n-1) / qchisq(alpha/(2*k), n-1))
    
    IC <- rbind(IC,
                ICmean,
                ICvar)
}
Ps
IC


##############
### K-MEANS 
#############

result.k <- kmeans(train_s, centers=4) # Centers: fixed number of clusters

names(result.k)

result.k$cluster      # labels of clusters
result.k$centers      # centers of the clusters
result.k$totss        # tot. sum of squares
result.k$withinss     # sum of squares within clusters
result.k$tot.withinss # sum(sum of squares nei cluster)
result.k$betweenss    # sum of squares between clusters
result.k$size         # dimention of the clusters

x11()
plot(data, col = result.k$cluster+1)

# in 3d (not useful for exam)
plot3d(Q, size = 3, col = result.k$cluster + 1, aspect = F)
points3d(result.k$centers, size = 10)

### How to choose k:
### 1) evaluate the variability between the groups with respect to 
###   the variability withing the groups
### 2) evaluate the result of hierarchical clustering (not recommended,
###    quite computationally expensive)

# method 1)
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

#method 2 expensive












