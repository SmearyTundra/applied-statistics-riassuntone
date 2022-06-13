diamonds <- read.table("diamonds.txt", header = T)
data <- scale(diamonds)
diss <- dist(data)
w <- hclust(diss, method = "ward.D2")
x11()
plot(w, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
#we have stable clustering structures with two or three clusters
#let's look at the induced partitions:
clust2 <- cutree(w, k = 2)
clust3 <- cutree(w, k = 3)

x11()
layout(cbind(1,2))
plot(data)
points(data[which(clust2 == 1),], col = 'red')
points(data[which(clust2 == 2),], col = 'blue')
plot(data)
points(data[which(clust3 == 1),], col = 'red')
points(data[which(clust3 == 2),], col = 'blue')
points(data[which(clust3 == 3),], col = 'green')

#by graphical inspection it seems clear how the two-clusters structure is the most appropriate
n1 <- length(which(clust2 == 1))
n2 <- length(which(clust2 == 2))
prob <- c(n1,n2)/(n1+n2)
prob

#we check if we have gaussianity of the data
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
mcshapiro.test(diamonds[which(clust2==1),])
mcshapiro.test(diamonds[which(clust2==2),])

#so at level 5% (one at the time) we cannot reject gaussianity for the two populations
#we use Bonferroni corrections on the non-scaled variables, to get more interpretability:
mean1 <- sapply(diamonds[which(clust2==1),], mean)
mean2 <- sapply(diamonds[which(clust2==2),], mean)
smean <- mean1-mean2
S1 <- cov(diamonds[which(clust2==1),])
S2 <- cov(diamonds[which(clust2==2),])
S1
S2
#same covariance structure
Sp <- ((n1 - 1)*S1 + (n2 - 1)*S2)/(n1+n2-2)

alphaB <- 0.05/2
band <- qt(1-alphaB/2, n1+n2-2)*sqrt(diag(Sp)*(1/n1 + 1/n2))
CI <- cbind(inf = smean - band,
            center = smean,
            sup = smean + band)
rownames(CI) <- c("Diam 1-2", "Car 1-2")
CI
#so the first (more populated cluster) is made of diamonds which are smaller and with less carates