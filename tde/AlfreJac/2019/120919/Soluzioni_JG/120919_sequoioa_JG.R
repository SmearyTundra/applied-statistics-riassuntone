sequoia <- read.table("sequoia.txt", header = T)
diss <- dist(sequoia)#euclidian by default
w <- hclust(diss, method = "ward.D2")
x11()
plot(w, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')

clust_2 <- cutree(w, k = 2)
clust_3 <- cutree(w, k = 3)
clust_4 <- cutree(w, k = 4)
clust <- cutree(w, k = 5)

x11()
layout(cbind(1,2,3,4))
plot(sequoia)
points(sequoia[which(clust_2==1),],col = 'red')
points(sequoia[which(clust_2==2),],col = 'blue')
plot(sequoia)
points(sequoia[which(clust_3==1),],col = 'red')
points(sequoia[which(clust_3==2),],col = 'blue')
points(sequoia[which(clust_3==3),],col = 'green')
plot(sequoia)
points(sequoia[which(clust_4==1),],col = 'red')
points(sequoia[which(clust_4==2),],col = 'blue')
points(sequoia[which(clust_4==3),],col = 'green')
points(sequoia[which(clust_4==4),],col = 'black')
plot(sequoia)
points(sequoia[which(clust==1),],col = 'red')
points(sequoia[which(clust==2),],col = 'blue')
points(sequoia[which(clust==3),],col = 'green')
points(sequoia[which(clust==4),],col = 'black')
points(sequoia[which(clust==5),],col = 'pink')

sequoia_1 <- sequoia[which(clust==2),]
sequoia_2 <- sequoia[which(clust==3),]
sequoia_3 <- sequoia[which(clust==1),]
sequoia_4 <- sequoia[which(clust==4),]
sequoia_5 <- sequoia[which(clust==5),]

centr_1 <- colMeans(sequoia_1)
centr_2 <- colMeans(sequoia_2)
centr_3 <- colMeans(sequoia_3)
centr_4 <- colMeans(sequoia_4)
centr_5 <- colMeans(sequoia_5)

n1 <- dim(sequoia_1)[1]
n2 <- dim(sequoia_2)[1]
n3 <- dim(sequoia_3)[1]
n4 <- dim(sequoia_4)[1]
n5 <- dim(sequoia_5)[1]



load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")

shapiro.test(sequoia_1$diameter)
shapiro.test(sequoia_2$diameter)
shapiro.test(sequoia_3$diameter)
shapiro.test(sequoia_4$diameter)
shapiro.test(sequoia_5$diameter)

alpha_B <- 0.1/(5*2)
ns <- c(n1,n2,n3,n4,n5)
q.tstd <- qt(1-alpha_B/2, df = ns-rep(1,5))
var <- (c(sd(sequoia_1$diameter),sd(sequoia_2$diameter), sd(sequoia_3$diameter), sd(sequoia_4$diameter), sd(sequoia_5$diameter) ))^2
means <- c(centr_1[2],centr_2[2],centr_3[2],centr_4[2],centr_5[2])
CI_mean <- cbind(means - q.tstd*sqrt(var/ns),
                 means,
                 means + q.tstd*sqrt(var/ns))
q.lowchi <- qchisq(alpha_B/2, df = ns-rep(1,5))
q.hichi <- qchisq(1-alpha_B/2, df = ns-rep(1,5))
CI_var <- cbind((ns-rep(1,5))*var/q.hichi,
                var,
                (ns-rep(1,5))*var/q.lowchi)
CI_mean
CI_var
