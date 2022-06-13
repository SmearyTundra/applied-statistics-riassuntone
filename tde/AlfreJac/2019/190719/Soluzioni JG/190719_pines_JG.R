pines <- read.table("pines.txt", header = T)
alpha <- 0.01
n <- dim(pines)[1]
p <- dim(pines)[2]


mu_0 <- c(14.2350,42.4520)

smean <- sapply(pines, mean)
S <- cov(pines)

#T_stat = (smean - mu_0)'(S/n)^-1(smean-mu_0) ~ (n-1)*p/(n-p) F(p, n - p)
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha, p, n-p)
T_stat <- n*t(smean - mu_0)%*%solve(S)%*%(smean-mu_0)

#H0: mean = mu_0 vs H1
Rej <- (T_stat > cfr.fisher)
Rej
p_val <- 1 - pf(T_stat/((n-1)*p/(n-p)), p, n - p)
p_val

load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
mcshapiro.test(pines)
#we can assume gaussianity, since it's gaussian in all its linear combinations


#This is a prediction ellipse, we should consider the variability of our data AROUND the mean
#analytical expression:
#   {x e R2 s.t. (x - mu_0)'(S)^-1(x - mu_0) <= radius^2} #if we want to assume that the true mean is mu_0 -> we have evidence supporting that
#
# (x - mean)'(SIGMA)^-1(x - mean) ~ Chi(p)
#ASYMPTOTICALLY => S -> SIGMA
#so radius = sqrt(quantile of a Chi)

radius <- sqrt(qchisq(1-alpha, p))

center <- mu_0
eigen(S)
dir1 <-  eigen(S)$vectors[,1]
dir2 <-  eigen(S)$vectors[,2]
length1 <- sqrt(eigen(S)$values[1])*radius
length2 <- sqrt(eigen(S)$values[2])*radius

x11()
plot(pines[,1],pines[,2], col='blue',pch=19,xlab='long',ylab='lat',asp=1)
points(mu_0[1],mu_0[2], pch=19,col = 'red')
points(smean[1],smean[2], pch=19,col = 'green')
ellipse(center=mu_0, shape=cbind(S), radius=radius, col = 'pink1')


#if we do not want to assume mu_0 to be the true center:
# (x - smean) ~ N(0, SIGMA + 1/n*SIGMA)
# (n/(n+1))*(x - smean)'(SIGMA)(x-smean) ~ Chi(2)
#asymptotically SIGMA =~ S

radius_new <- sqrt(((n+1)/n)*qchisq(1-alpha, p))

center_new <- smean
eigen(S)
dir1 <-  eigen(S)$vectors[,1]
dir2 <-  eigen(S)$vectors[,2]
length1_new <- sqrt(eigen(S)$values[1])*radius_new
length2_new <- sqrt(eigen(S)$values[2])*radius_new

x11()
plot(pines[,1],pines[,2], col='blue',pch=19,xlab='long',ylab='lat',asp=1)
points(mu_0[1],mu_0[2], pch=19,col = 'red')
points(smean[1],smean[2], pch=19,col = 'green')
ellipse(center=smean, shape=cbind(S), radius=radius_new, col = 'pink1')

#confronting with the previous one
ellipse(center=mu_0, shape=cbind(S), radius=radius, col = 'red')
#quite similar