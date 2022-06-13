Sailing <- read.table("Sailing.txt",header = T)
attach(Sailing)
#If the two populations have same covariance structure we can use LDA
levels(as.factor(type))
l <- lda(type ~ water + sailing.time, prior = c(0.2,0.8))
l
#if the two populations might have different covariance structure we MUST use QDA
q <- qda(type ~ water + sailing.time, prior = c(0.2,0.8))
q

#plot LDA
w <- seq(min(water),max(water), length = 200)
s <- seq(min(sailing.time),max(sailing.time), length = 200)
ws <- expand.grid(water=w, sailing.time=s)
zl  <- predict(l, ws)$post

zl1 <- zl[,1] - zl[,2]
#zl1 <- zl[,1] - pmax(zl[,2], zl[,3])
#zl2 <- zl[,2] - pmax(zl[,1], zl[,3])

x11()
plot(water, sailing.time, main = "LDA")
points(water[which(type == "seadog")], sailing.time[which(type == "seadog")], col = 'blue')
points(water[-which(type == "seadog")], sailing.time[-which(type == "seadog")], col = 'red')
legend("topright", legend=levels(as.factor(type)), fill=c('blue','red'), cex=.7)
points(l$means, pch=4,col=c('blue','red') , lwd=2, cex=1.5)
contour(w, s, matrix(zl1, 200), levels=0, drawlabels=F, add=T)  
#contour(w, s, matrix(zl2, 200), levels=0, drawlabels=F, add=T)  

zq  <- predict(q, ws)$post
zq1 <- zq[,1] - zq[,2]

x11()
plot(water, sailing.time, main = "QDA")
points(water[which(type == "seadog")], sailing.time[which(type == "seadog")], col = 'blue')
points(water[-which(type == "seadog")], sailing.time[-which(type == "seadog")], col = 'red')
legend("topright", legend=levels(as.factor(type)), fill=c('blue','red'), cex=.7)
points(l$means, pch=4,col=c('blue','red') , lwd=2, cex=1.5)
contour(w, s, matrix(zq1, 200), levels=0, drawlabels=F, add=T)  


post_l <- predict(l, Sailing[1,1:2])$post
post_q <- predict(q, Sailing[1,1:2])$post


length(which(type=="seadog"))/dim(Sailing)[1]
#the prior probabilities are not the empirical ones, obtained from the dataset
#APER does not provide an estimate of AER
#we need to adapt it

lass <- predict(l, Sailing[,1:2])
qass <- predict(q, Sailing[,1:2])

prior <- c(0.2, 0.8)
G <- 2
miscl <- table(class.true=type, class.assigned=lass$class)
 APERl <- 0
 for(g in 1:G)
   APERl <- APERl + sum(miscl[g,-g])/sum(miscl[g,]) * prior[g]  
APERl
 
miscq <- table(class.true=type, class.assigned=qass$class)
APERq <- 0
for(g in 1:G)
  APERq <- APERq + sum(miscq[g,-g])/sum(miscq[g,]) * prior[g]  
APERq

#they yield pretty much the same performance


z_new <- data.frame(water = 35, sailing.time = 168)
l_pred <- predict(l, z_new)$post
q_pred <- predict(q, z_new)$post

S1 <- cov(Sailing[which(type=="seadog"), 1:2])
S2 <- cov(Sailing[-which(type=="seadog"), 1:2])

var.test(sailing.time ~ type)
bartlett.test(sailing.time ~ type)
bartlett.test(Sailing[,1:2],type)

#we cannot assume homoschedasticity -> better to use the prediction yielded by QDA
#The new observation should be classified as a ship full of peasants
