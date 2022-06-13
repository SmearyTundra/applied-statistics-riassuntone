tourists <- read.table("tourists.txt", header = T)
data <- tourists[,3:10]
n <- dim(tourists)[1]
x11()
barplot(sapply(data,sd)^2, las=2, main='Original variables', ylab='Variances')
#variances are not comparable, we scale the variables
data <- data.frame(scale(data))
pc <- princomp(data, scores = T)
summary(pc)
x11()
barplot(sapply(data,sd)^2, las=2, main='Original variables', ylab='Variances')

pc$loadings
x11()
barplot(pc$loadings[,1], main = "1st")
library(mvtnorm)
abline(h=0, v=0, lty=2, col='grey')
#all scores are almost equal -> basically a mean

x11()
barplot(pc$loadings[,2], main = "2nd")
abline(h=0, v=0, lty=2, col='grey')
#it's a delta between "luxury" and "cheap" sojourns
x11()
barplot(pc$loadings[,3], main = "3rd")
abline(h=0, v=0, lty=2, col='grey')
sum(pc$loadings[,3])# -> 0.06  it's actually balanced'
#it's basically a delta between bed and breakfast vs everthing else

point_a_loads <- pc$loadings[,1:3]
#the percentage of variance explained is found on "cumulative variance" on summary

proj <- pc$scores[,1:2]
x11()
plot(proj)
abline(h=0, v=0, lty=2, col='grey')

#1st quadrant -> high influx of tourists, mainly on luxury sojourns
#2nd/3rd quadrants -> low influx of tourists,not much economically differentiated, but in the 2nd quadrant they should prefer more luxurious sojourn than in the 3rd
#4th -> low influx of tourists, mainly on cheap sojourns

#the shape suggests that as the flux increases, the distinction between luxury and cheap locations increases as well
new_proj <- cbind(tourists[,1:2], proj)
colo <- new_proj$Month
k = 1
for(month in unique(tourists$Month)){
  colo[which(colo==month)] <- rep(rainbow(12)[k], length(which(colo==month)))
  k = k + 1
}
x11()
plot(new_proj[,3:4], col = colo)
legend('topright',levels(as.factor(tourists[,1])),fill=rainbow(12),bty='n')
#no particular patterns

colo <-  new_proj$Region.of.origin
k = 1
for(region in unique(tourists$Region.of.origin)){
  colo[which(colo==region)] <- rep(rainbow(length(unique(tourists$Region.of.origin)))[k], length(which(colo==region)))
  k = k + 1
}
x11()
plot(new_proj[,3:4], col = colo)
legend('topright',levels(as.factor(tourists[,2])),fill=rainbow(length(unique(tourists$Region.of.origin))),bty='n')
#maybe some slight dependencies on the region of provenience

x11()
boxplot(new_proj[,3]~new_proj$Region.of.origin)#Indeed there are differences in the groupings, completely explainable by the different number of residents in each region
boxplot(new_proj[,4]~new_proj$Region.of.origin)

boxplot(new_proj[,3]~new_proj$Month)
boxplot(new_proj[,4]~new_proj$Month)#as we suspected there is no evident relationship wrt to month
summary(pc)
#dimensionality reduction -> 1st PC, just by itself captures 89% of the variability
#for a more precise picture we can consider the 2nd PC as well, still interpretable and quite significant; it's useless to go further than that anyway