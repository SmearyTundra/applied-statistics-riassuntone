sunchair <- read.table("sunchair.txt", header = T)
attach(sunchair)
n <- dim(sunchair)[1]
g <- dim(sunchair)[2]
n_tot <- g*n
#we need to test for gaussianity amongst the four groups
shapiro.test(Mar.May)
shapiro.test(Jun.July)
shapiro.test(Aug.Oct)
shapiro.test(Nov.Feb)
#We can assume gaussianity
#we need to test for the same variance
#graphically:
x11()
boxplot(sunchair)
#not ideal, but we can work on that

#surely there is a better method, but the best that I can come up with is:
data <- data.frame(c(Mar.May, Jun.July,Aug.Oct,Nov.Feb))
labels <- as.factor(c(rep("Mar.May", n), rep("Jun.July", n), rep("Aug.Oct", n), rep("Nov.Feb", n)))
data<- cbind(prices = data, labels)
names(data) <- c("price", "period")
an <- aov(data$price ~ data$period)
summary(an)

#there is a significant variation of prices during the year

SSres <- sum((an$residuals)^2)
#Test statistic: (delta_sample_mean - delta_mean)/(sqrt(SSres/(n_tot-g)*(1/n_i + 1/n_j))) ~ t(n_tot - g)
alpha_bonf <- 0.05/4
rf.T <- qt(1-alpha_bonf/2, n_tot - g)
#since all groups have the same number of observations we can simplify:
means <- colMeans(sunchair)
centers <- matrix(c(means[2] - means[1],means[3] - means[2],means[4] - means[3], means[1] - means[4]))
increment <- rf.T*sqrt(SSres/(n_tot-g)*(2/n))
CI <- cbind(centers - rep(increment, g), centers, centers + rep(increment,g))
colnames(CI)<-c("inf", "point", "sup")
rownames(CI)<-c("Jun.July - Mar.May", "Aug.Oct - Jun.July", "Nov.Feb - Aug.Oct", "Mar.May - Nov.Feb")
CI
#so the best period to buy sunchairs is November->February








#in case dynamics is not interpreted as "show the differences of the meanover time", but as "show how the mean changes over time", which is actually a better interpretation of the exercise:

#Test statistic: (sample_mean - mean)/(sqrt(S/n)) ~ t(n-1)
alpha_bonf <- 0.05/4
S1 <- var(Mar.May)
S2 <- var(Jun.July)
S3 <- var(Aug.Oct)
S4 <- var(Nov.Feb)
rf.T <- qt(1-alpha_bonf/2, n - 1)
means <- matrix(colMeans(sunchair))
increment <- matrix(rf.T*sqrt(c(S1,S2,S3,S4)/n))
CI <- cbind(means - increment, means, means + increment)
colnames(CI)<-c("inf", "point", "sup")
rownames(CI)<-c("Mar.May", "Jun.July", "Aug.Oct", "Nov.Feb")
CI

#same interpretation -> Nov.Feb is the best period to buy a sunchair