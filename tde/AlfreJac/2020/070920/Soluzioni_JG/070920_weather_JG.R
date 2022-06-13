weather <- read.table("weather.txt", header = T)
x11()
barplot(sapply(weather,sd)^2, las=2, main='Original variables', ylab='Variances')
#not comparable sample variances

pc_nstd <- princomp(weather, scores = T)
summary(pc_nstd)
x11()
barplot(pc_nstd$loadings[,1])
x11()
barplot(pc_nstd$loadings[,2])

#indeed we get a structure very similar to the original variances of our dataset

#we use standardized variables

data <- scale(weather,scale = T, center = T)
pc <- princomp(data, scores = T)
summary(pc)

x11()
barplot(diag(cov(data)))

x11()
barplot(pc$loadings[,1], ylim = c(-1,1), main = "PC1")
#basically an average on the indicatores of temperature -> the higher the score the higher the temperature
x11()
barplot(pc$loadings[,2], ylim = c(-1,1), main = "PC2")
#basically an average on wind indicators -> the higher the score the slower the winds
x11()
barplot(pc$loadings[,3], ylim = c(-1,1), main = "PC3")
#difference between humidity and visibility -> the higher the more humid and the lower the visibility is

x11()
plot(pc$scores[,1], pc$scores[,2], main = "PC1 & PC2 plane", xlab = "Scores-1", ylab = "Scores-2")
abline(h = 0, v = 0, lty = 2, col = "grey")
colournames <- rainbow(n = dim(weather)[1],start = 0.45, end = 0.8)
x11()
plot(pc$scores[,1], pc$scores[,2], main = "PC1 & PC2 plane", xlab = "Scores-1", ylab = "Scores-2", col = colournames)

#scree plot
x11()
layout(cbind(c(1,3), c(2,2)))
barplot(pc$sd)
plot(cumsum(pc$sdev)/sum(pc$sd), type = 'l', ylim = c(0,1))
abline(h = 0.8, lty = 2, col = "blue")
barplot(diag(cov(data)))

#we can just consider the first tree principal components, since they account for 94% of the variability
summary(pc)
#we have the proportion of variability on the summary


#projections : scalar prod: z'PC1, z'PC2, z'PC3 -> of the SCALED unit!!!
z <- c(30, 23, 36, 22, 65, 19, 5, 15)
scale_z <- (z - colMeans(weather))/sapply(weather,sd)
proj <- sum(scale_z*pc$loadings[,1])*pc$loadings[,1] + sum(scale_z*pc$loadings[,2])*pc$loadings[,2] + sum(scale_z*pc$loadings[,3])*pc$loadings[,3]
coeff <- c(sum(scale_z*pc$loadings[,1]), sum(scale_z*pc$loadings[,2]), sum(scale_z*pc$loadings[,3]))