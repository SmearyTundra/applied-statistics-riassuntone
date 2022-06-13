####
# PROBLEMA 1
####

data <- read.table('kimono.txt')
head(data)
dim(data)

data$city <- factor(data$city)
data$type <- factor(data$type)

city <- levels(data$city)
type <- levels(data$type)

i1 <- which(data$type == 'hand-made' & data$city== 'Kyoto')
i2 <- which(data$type == 'ready-to-use' & data$city== 'Kyoto')
i3 <- which(data$type == 'hand-made' & data$city== 'Tokyo')
i4 <- which(data$type == 'ready-to-use' & data$city== 'Tokyo')

#assumptions
shapiro.test(data[i1,1])
shapiro.test(data[i2,1])
shapiro.test(data[i3,1])
shapiro.test(data[i4,1])

bartlett.test((list(data[i1,1],
                    data[i2,1],
                    data[i3,1],
                    data[i4,1])))

fit.aov <- aov(value ~ city + type + city:type, data)
summary.aov(fit.aov)

fit.aov <- aov(value ~ city + type, data)
summary.aov(fit.aov)

fit.aov <- aov(value ~  type, data)
summary.aov(fit.aov)
attach(data)
Mediag  <- tapply(value, type, mean) #menas in each group 
SSres <- sum(residuals(fit.aov)^2)
n <- dim(data)[1]
g <- 2
S <- SSres/(n-g)
n1 <- length(i1)+ length(i3)
n2 <- length(i2)+ length(i4)
paste(levels(type)[1],"-",levels(type)[2])
alpha <- 0.05
as.numeric(c(Mediag[1]-Mediag[2] - qt(1-alpha/(2), n-g) * sqrt( S * ( 1/n1 + 1/n2 )),
             Mediag[1]-Mediag[2] + qt(1-alpha/(2), n-g) * sqrt( S * ( 1/n1 + 1/n2 ))))


######
# PROBLEMA 2
#####
data <- read.table('bento.txt')
head(data)
dim(data)
D <- data[,1:4] -data [,5:8]
head(D)
load("~/Desktop/IV anno/Applied statistics/mcshapiro.test.RData")
mcshapiro.test(D)

n <- dim(D)[1]  
p <- dim(D)[2]  

D.mean   <- sapply(D,mean)
D.cov    <- cov(D)
D.invcov <- solve(D.cov)

delta.0 <- c(0,0,0,0)

D.T2 <- n * (D.mean-delta.0) %*% D.invcov %*% (D.mean-delta.0) #test statistic(Malanobius distance)
D.T2

# we compute the p-value
P <- 1-pf(D.T2*(n-p)/(p*(n-1)), p, n-p)
P
alpha <- 0.05
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p) #radius of the ellipse
cfr.fisher


D <- -D
IC.T2.1 <- c( D.mean[1]-sqrt(cfr.fisher*D.cov[1,1]/n) , D.mean[1], D.mean[1]+sqrt(cfr.fisher*D.cov[1,1]/n) )
IC.T2.2  <- c( D.mean[2]-sqrt(cfr.fisher*D.cov[2,2]/n) , D.mean[2], D.mean[2]+sqrt(cfr.fisher*D.cov[2,2]/n) )
IC.T2.3 <- c( D.mean[3]-sqrt(cfr.fisher*D.cov[3,3]/n) , D.mean[3], D.mean[3]+sqrt(cfr.fisher*D.cov[3,3]/n) )
IC.T2.4  <- c( D.mean[4]-sqrt(cfr.fisher*D.cov[4,4]/n) , D.mean[4], D.mean[4]+sqrt(cfr.fisher*D.cov[4,4]/n) )

T2 <- rbind(IC.T2.1, IC.T2.2, IC.T2.3 , IC.T2.4)
dimnames(T2)[[2]] <- c('inf','center','sup')
T2

#####
# PROBLEMA 3
#####

data <- read.table('geisha.txt')
head(data)
dim(data)
x11()
plot(data)
# we have the four features but non the lablel -> we can see that there are two separated grups 
euc <- dist(data, method='euclidean')
es <- hclust(euc, method='average')
x11()
plot(es, main='euclidean-avarage', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
cluster.es <- cutree(es, k=2)
x11()
plot(data, col=ifelse(cluster.es==1,'red','blue'), pch=19)

coph.es <- cophenetic(es)
es <- cor(euc, coph.es)
es
C <- colMeans(data[cluster.es == 1,])
# compute the centroid
C <- colMeans(data[cluster.es == 2,])
C
dim(data)[1]

i1 <- which(cluster.es == 1)
i2 <- which(cluster.es == 2)

load("~/Desktop/IV anno/Applied statistics/mcshapiro.test.RData")
mcshapiro.test(data[i1,])
mcshapiro.test(data[i2,])
cov(data[i1,])
cov(data[i2,])
cluster.es
Mean1 <- sapply(data[i1,],mean)
Mean2 <- sapply(data[i2,],mean)
S1 <- cov(data[i1,])
n1 <- length(i1)
S2 <- cov(data[i2,])
n2 <- length(i2)
alpha <- 0.05
k<- 4
pivotal_diff <- qnorm(1-alpha/k)
pivotal_mean <- qt(1-alpha/k , n1-1)
IC.1 <- c(Mean1[1] - Mean2[1] - pivotal_diff* sqrt((S1/n1+ S2/n2)[1,1]),
          Mean1[1] - Mean2[1] + pivotal_diff* sqrt((S1/n1+ S2/n2)[1,1]))
IC.2 <- c(Mean1[2] - Mean2[2] - pivotal_diff* sqrt((S1/n1+ S2/n2)[2,2]),
          Mean1[2] - Mean2[2] + pivotal_diff* sqrt((S1/n1+ S2/n2)[2,2]))
IC.4 <- c(Mean1[2] - pivotal_mean* sqrt((S1/n1)[2,2]),
          Mean1[2] + pivotal_mean* sqrt((S1/n1)[2,2]))
IC.3 <- c(Mean1[1] - pivotal_mean* sqrt((S1/n1)[1,1]),
            Mean1[1] + pivotal_mean* sqrt((S1/n1)[1,1]))
IC <- rbind(IC.1, IC.2, IC.3, IC.4)
IC

######
#PROBLEM 4
#####
data <- read.table('garden.txt')
head(data)
dim(data)

fm <- lm(extension ~ .,data)

summary(fm) 
x11()
par(mfrow=c(2,2))
plot(fm)
shapiro.test(residuals(fm))

fm$coefficients
library(MASS)
library(car)
library(rgl)

linearHypothesis(fm, rbind(c(0,0,1,0,0), c(0,0,0,1,0)), c(0,0))
linearHypothesis(fm, rbind(c(0,1,0,0,0), c(0,0,0,0,1)), c(0,0))

vif(fm)

cov(data[,1:4])
data.pc <- princomp(data[,1:4], scores=TRUE)
# primcomp by default center the data 
summary(data.pc)

par(mar = c(1,4,0,2), mfrow = c(2,1))
for(i in 1:2)barplot(data.pc$load[,i], ylim = c(-1, 1))

dev.off()


sp1.pc <- data.pc$scores[,1]
sp2.pc <- data.pc$scores[,2]
extention <- data [,5]
fm.pc <- lm(extention ~ sp1.pc + sp2.pc)

summary(fm.pc) 
attach(data)
m1 <- mean(carps)
m2 <- mean(maple)
m3 <- mean(cherry)
m4 <- mean(stones)
beta0 <- coefficients(fm.pc)[1] - 
  coefficients(fm.pc)[2]*data.pc$load[1,1]*m1 - 
  coefficients(fm.pc)[3]*data.pc$load[1,2]*m1 - 
  coefficients(fm.pc)[2]*data.pc$load[2,1]*m2 - 
  coefficients(fm.pc)[3]*data.pc$load[2,2]*m2 -
  coefficients(fm.pc)[2]*data.pc$load[3,1]*m3 - 
  coefficients(fm.pc)[3]*data.pc$load[3,2]*m3 -
  coefficients(fm.pc)[2]*data.pc$load[4,1]*m4 - 
  coefficients(fm.pc)[3]*data.pc$load[4,2]*m4
beta1 <- coefficients(fm.pc)[2]*data.pc$load[1,1] + 
  coefficients(fm.pc)[3]*data.pc$load[1,2]  
beta2 <- coefficients(fm.pc)[2]*data.pc$load[2,1] + 
  coefficients(fm.pc)[3]*data.pc$load[2,2] 
beta3 <- coefficients(fm.pc)[2]*data.pc$load[3,1] + 
  coefficients(fm.pc)[3]*data.pc$load[3,2]  
beta4 <- coefficients(fm.pc)[2]*data.pc$load[4,1] + 
  coefficients(fm.pc)[3]*data.pc$load[4,2] 

