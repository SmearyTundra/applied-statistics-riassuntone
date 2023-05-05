###########################################################################
##### PCA PRINCIPAL COMPONENT ANALYSIS
###########################################################################

library(car)

##### Pca + Projection (No standard) ####

weather <- read.table('weather.txt', header = TRUE)
head(weather)

n <- dim(weather)[1]
p <- dim(weather)[2]

# Boxplot
par(mar=rep(8,4))
boxplot(weather, las=2, col='gold')

# Note: PCA is not about the mean, it is about the variability
par(mar=rep(8,4))
boxplot(scale(x=weather,center = T, scale=F), las=2, col='gold')

pc.weather <- princomp(weather, scores=T)
pc.weather
summary(pc.weather)

# loadings (recall: coefficients of the linear combination of the original 
#           variables that defines each principal component)

load.wea <- pc.weather$loadings
load.wea

load.wea[,1:8]

biplot(pc.weather)

# graphical representation of the loadings of the first 3 principal components
par(mfcol = c(3,1))
for(i in 1:3) barplot(load.wea[,i], ylim = c(-1, 1), main=paste("PC",i))

#Explained variance (screeplot)
layout(matrix(c(2,3,1,3),2,byrow=T))
plot(pc.weather, las=2, main='Principal components', ylim=c(0,4.5e2))
barplot(sapply(weather,sd)^2, las=2, main='Original Variables', ylim=c(0,4.5e2), ylab='Variances')
plot(cumsum(pc.weather$sd^2)/sum(pc.weather$sd^2), type='b', axes=F, xlab='number of components', 
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(weather),labels=1:ncol(weather),las=2)

#Scatterplot
# scores
scores.wea <- pc.weather$scores
scores.wea

plot(scores.wea[,1:2])
abline(h=0, v=0, lty=2, col='grey', pch = 16)

day <- seq(1,213,1)
plot(scores.wea[,1:2], col=ifelse(day < 90, 'red', ifelse (day < 160, 'green', 'black') ), pch = 16)
abline(h=0, v=0, lty=2)

#projection
new <- data.frame( MeanTemp = 30, MinTemp = 23, MaxTemp = 36, DewPoint = 22, Humidity = 65, Visibility = 19, MeanWind = 5, MaxWind = 15)
scores.new <- t(pc.weather$loadings)%*%as.numeric(new-colMeans(weather))
scores.new

plot(scores.wea[,1:2])
abline(h=0, v=0, lty=2, col='grey')
points(scores.new[1],scores.new[2], col = 'red', pch = 19)






# projection with rescaled dataset
x.mean <- colMeans(data)
x.var <- sapply(data, FUN = sd)

new <- c(
    Energy_kcal=400,
    Protein_g=9,
    Fat_g=5,
    Carb_g=100,
    Sugar_g=30,
    Fiber_g=4
)

new <- (new - x.mean) / x.var

pc.proj <- new %*% pc$loadings[,1:2]
pc.proj

par(mfrow=c(1,1))
plot(pc$scores)
points(pc.proj, col=2, pch=16, cex=2)