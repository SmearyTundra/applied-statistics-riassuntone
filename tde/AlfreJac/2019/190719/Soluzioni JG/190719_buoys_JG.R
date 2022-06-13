buoys <- read.table("buoys.txt", header = T)
data <- buoys[,1:2]
diss <- dist(data)#euclidian by default

hcw <- hclust(diss, method = "ward.D2")
x11()
plot(hcw)
#we can either use two or three clusters, let's evaluate graphically
x11()
plot(data)
#we may have found the tree clusters
clust <- cutree(hcw, k = 3)
points(data[which(clust == 1),], col = 'red')
points(data[which(clust == 2),], col = 'blue')
points(data[which(clust == 3),], col = 'green')
#yep we work with three clusters
clust <- as.factor(clust)

data1 <- data[which(clust == 1),]
data2 <- data[which(clust == 2),]
data3 <- data[which(clust == 3),]

center1 <- sapply(data1, mean)
center2 <- sapply(data2, mean)
center3 <- sapply(data3, mean)

n1 <- dim(data1)[1]
n2 <- dim(data1)[2]
n3 <- dim(data1)[3]



#we perform a test anova, hence considering x = mu + tau_i + eps; i = 1,2,3; eps ~ iid N(0, sigma^2)
#we can already assume independence on the units
#we need to verify that the variability amongst groups is the same and that the residuals have a gaussian behaviour
#the test is H0: tau_1 = tau_2 = tau_3 = 0 vs H1

an <- aov(buoys$DO ~ clust)
summary(an)

shapiro.test(an$residuals)
shapiro.test(an$residuals[which(clust==1)])
shapiro.test(an$residuals[which(clust==2)])
shapiro.test(an$residuals[which(clust==3)])

x11()
qqnorm(an$residuals)
qqline(an$residuals)
x11()
plot(an$fitted.values, an$residuals)

var.test(buoys$DO[which(clust==1)],buoys$DO[which(clust==2)])
var.test(buoys$DO[which(clust==1)],buoys$DO[which(clust==3)])
var.test(buoys$DO[which(clust==3)],buoys$DO[which(clust==2)])

#the assumptions of the model have been verified


#there is no statistical evidence to support the claim of a difference in the mean of DO
#indeed this can also be clearly seen in
x11()
boxplot(buoys$DO ~ clust)
