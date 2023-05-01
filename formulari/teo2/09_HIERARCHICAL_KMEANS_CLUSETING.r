### TOPICS:
### Hierarchical clustering
### K-means clustering
### Exercises
### Multidimensional Scaling
### Plotting region containing 99&% of data

library(mvtnorm)
library(rgl)
library(car)





############### HIERARCHICAL CLUSTERING ###################
###########-------------------------------------------------

n <- nrow(data) # number of samples

# d = distance
# l = linkage

#### COMPUTING DISTANCE
###-------------------
# distance matrix
dist.matrix <- dist(data, method='euclidean-canberra-manhattan')

# plot dist matrix
image(1:n, 1:n, as.matrix(dist.matrix), main = "metrics: NAME OF METRIC", asp = 1, xlab = "i", ylab = "j")

# clustering
cluster.dl <- hclust(dist.matrix, method = "single-average-complete-ward.D2")
names(cluster.dl)
# merge  # order of aggregation of statistical units / clusters
# height # distance at which we have aggregations
# order  # ordering that allows to avoid intersections in the dendrogram

#### DENDROGRAM
###-------------------
plot(cluster.dl, main = "distance-linkage", hang = -0.1, xlab = "", labels = F, cex = 0.6, sub = "")
# with k cluster
rect.hclust(cluster.dl, k = 2-3-4)

#### CUTTING
###-------------------
cluster.dl.cut <- cutree(cluster.dl, k = 2)

#### COPH COEFFICIENTS
###-------------------
# cophenetic matrix
coph.matrix <- cophenetic(cluster.dl)
# cophenetic coefficient
coph.coeff <- cor(dist.matrix, coph.matrix)
coph.coeff

#### NUMEROSITY OF CLUSTERS
###-------------------
table(cluster.dl.cut)

#### CENTROIDS
###-------------------
centroids <- apply(data, dim(data)[2], function (x) tapply (x, cluster.es, mean))
centroids

#### PLOT DISTINCTION
###-------------------
plot(data, col = cluster.dl.cut+1, pch = 19) # or col = cluster.dl (+1 for changing colors)









#### BONFERRONI CONFIDENCE FOR CLUSTERING (+ gaussianity)
###-------------------

# 1) univariate, just the mean/variable of one variable in different clusters
# MEDIA + VARIANZA
alpha <- 0.05
k <- 3 # number of CI
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



# 2) Bivariate, mean of two variables in different clusters (no variance)
alpha <- 0.05
k <- 6 # number of CI
g <- 3 # number of clusters
IC={}
Ps={}
for(i in 1:g){
    X <- d[cluster.es==i,1] # min
    Y <- d[cluster.es==i,2] # artists
    n <- length(X)
    Ps <- c(Ps,shapiro.test(X)$p)
    Ps <- c(Ps,shapiro.test(Y)$p)
    x.mean <- mean(X)
    y.mean <- mean(Y)
    x.cov <- var(X)
    y.cov <- var(Y)
    
    ICmeanMin <- c(inf    = x.mean - sqrt(x.cov/n) * qt(1 - alpha/(2*k), n-1),
                   center = x.mean,
                   sup    = x.mean + sqrt(x.cov/n) * qt(1 - alpha/(2*k), n-1))
    
    ICmeanArtist <- c(inf    = y.mean - sqrt(y.cov/n) * qt(1 - alpha/(2*k), n-1),
                      center = y.mean,
                      sup    = y.mean + sqrt(y.cov/n) * qt(1 - alpha/(2*k), n-1))
    
    IC <- rbind(IC,
                ICmeanMin,
                ICmeanArtist)
}
Ps # Assumptions: normalità verificata -> Ps
IC



















############### K-MEANS CLUSTERING ###################
###########-------------------------------------------------
## in automatic, command kmeans()
result.k <- kmeans(Q, centers=2) # Centers: fixed number of clusters
names(result.k)

result.k$cluster      # labels of clusters
result.k$centers      # centers of the clusters
result.k$totss        # tot. sum of squares
result.k$withinss     # sum of squares within clusters
result.k$tot.withinss # sum(sum of squares within cluster)
result.k$betweenss    # sum of squares between clusters
result.k$size         # dimension of the clusters

x11()
plot(Q, col = result.k$cluster+1)


#### CHOOSE K
### --------
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















#### Region containing 99% of data #####
###-------------------------------------
mcshapiro.test(p.pr) # if true

M <- sapply(p.pr, mean)
S <- cov(p.pr)
alpha <- 0.01
cfr.chisq <- qchisq(1-alpha,p)

# Characterize the ellipse:
# Axes directions:
eigen(S)$vectors
# Center:
M
# Radius of the ellipse:
r <- sqrt(cfr.chisq)
# Length of the semi-axes:
r*sqrt(eigen(S)$values)  

x11()
plot(p.pr, asp = 1, col='gold', pch=19, xlim=c(-10,50))
points(M[1], M[2], pch = 4, cex = 1.5, lwd = 2)
ellipse(center=M, shape=S, radius=sqrt(cfr.chisq), col = 'black', lty = 2, center.pch = 4)
points(stones[which(cl.ew==1),])