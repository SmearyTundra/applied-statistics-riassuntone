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
species.name <- d[,5]
d        <- d[,1:4]


#### COMPUTING DISTANCE
###-------------------
d.e <- dist(d, method='euclidean')
d.m <- dist(d, method='manhattan')
d.c <- dist(d, method='canberra')

x11()
par(mfrow=c(1,3))
image(1:150,1:150,as.matrix(d.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j' )
image(1:150,1:150,as.matrix(d.c), main='metrics: Canberra', asp=1, xlab='i', ylab='j' )
image(1:150,1:150,as.matrix(d.m), main='metrics: Manhattan', asp=1, xlab='i', ylab='j' )

# Unorder data
misc <- sample(n)
d <- d[misc,]
d.e <- dist(d, method='euclidean')
d.m <- dist(d, method='manhattan')
d.c <- dist(d, method='canberra')



#### CLUSTERING
###-------------------
d.es <- hclust(d.e, method='single')
d.ea <- hclust(d.e, method='average')
d.ec <- hclust(d.e, method='complete')
clustw <- hclust(d.e, method='ward.D2')
plot(clustw, hang=-0.1, labels=FALSE, main='ward', xlab='', sub='')


# if we want more detailed information on euclidean-complete
# clustering:
names(d.ec)
d.ec$merge  # order of aggregation of statistical units / clusters
d.ec$height # distance at which we have aggregations
d.ec$order  # ordering that allows to avoid intersections in the dendrogram




#### DENDROGRAMS
###-------------------
x11()
par(mfrow=c(1,3))
plot(d.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(d.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(d.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')

# Dendrograms 2 clusters etc
x11()
par(mfrow=c(1,3))
plot(d.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(d.es, k=2)
plot(d.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(d.ec, k=2)
plot(d.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(d.ea, k=2)
#### Look also at the scatter plot before deciding where to cut




#### CUTTING
###-------------------
# Fix k=2 clusters:
cluster.ec <- cutree(d.ec, k=2) # euclidean-complete:
cluster.es <- cutree(d.es, k=2) # euclidean-single
cluster.ea <- cutree(d.ea, k=2) # euclidean-average


# Interpretation?
# Interpret the clusters
table(label.true = species.name[misc], label.cluster = cluster.es)  # MISC SE HO RIMISCHIATO
table(label.true = species.name[misc], label.cluster = cluster.ec)
table(label.true = species.name[misc], label.cluster = cluster.ea)



#### COPH COEFFICIENTS
###-------------------
coph.es <- cophenetic(d.es)
coph.ec <- cophenetic(d.ec)
coph.ea <- cophenetic(d.ea)

# compute cophenetic coefficients
es <- cor(d.e, coph.es)
ec <- cor(d.e, coph.ec)
ea <- cor(d.e, coph.ea)

c("Eucl-Single"=es,"Eucl-Compl."=ec,"Eucl-Ave."=ea)




#### NUMEROSITY OF CLUSTERS
###-------------------
g <- 2
i1 <- which(cluster.ec==1)
i2 <- which(cluster.ec==2)
ng <- c(length(i1),length(i2)) 
ng 
N <- sum(ng)

# oppure
table(cluster.ec)





#### CENTROIDS
###-------------------
centroids <- apply(d, dim(d)[2], function (x) tapply (x, cluster.es, mean))
centroids




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
    X <- d[cluster.dl.cut==i,1] # 1 seleziona mean per prima var, 2 per la seconda
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




