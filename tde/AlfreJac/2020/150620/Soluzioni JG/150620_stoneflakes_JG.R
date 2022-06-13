stoneflakes <- read.table("stoneflakes.txt", header = T)
attach(stoneflakes)
diss <- dist(stoneflakes)#euclidian by default
wl <- hclust(diss, method = "ward.D2")
x11()
plot(wl, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
#by looking at the dendogram it is evident that the most stable clustering structure involves three clusters
#this intuition is confirmed by looking at the plot
x11()
plot(stoneflakes)

clust <- as.factor(cutree(wl, k = 3))
load("C:/Users/jacop/Desktop/università/da dare/Applied Statistics/AS lab/LAB_5/mcshapiro.test.RData")
#the model is:
#x = [Length, Width],
# x_i,j = tau_i + eps_i,j, i = 1,2,3; eps_i.j ~ N(0, SIGMA)

#we are testing
#H0: tau_1 = tau_2 = tau_3 vs H1

man <- manova(cbind(Length, Width) ~ clust)


summary(man)
summary.aov(man)
#there is evidence supporting the fact that at least two of the three different clusters have different means
#in order to say that we need to verify gaussianity of eps_i,j, since we've assume homoschedasticity
mcshapiro.test(man$residuals)
#we do not have evidence to reject the gaussian hypothesis, we can work with it

#we want confidence intervals for:
# Length_1 - Length_2
# Length_1 - Length_3
# Length_2 - Length_3
# Width_1 - Width_2
# Width_1 - Width_3
# Width_2 - Width_3
#so we want 6 simultaneous confidence intervals
n_1 <- length(which(clust==1))
n_2 <- length(which(clust==2))
n_3 <- length(which(clust==3))
n_1
n_2
n_3
#the groups do not have the same number of elements

#the general test statistic is:
# x_bar_i ~ N(mu_i, SIGMA/n_i)
# delta_bar_i,j ~ N(mu_i - mu_j, SIGMA(1/n_i + 1/n_j))
# delta_bar_i,j[k] ~ N(mu_i[k] - mu_j[k], SIGMA(1/n_i + 1/n_j)[k,k])

# (delta_bar_i,j[k] - truedelta)/sqrt(SIGMA(1/n_i + 1/n_j)[k,k]) ~ N(0, 1)
# S_pooled*(n_tot - g) ~ W(SIGMA, n_tot - g)
# k'S_pooledk*(n_tot - g)/SIGMA[k,k] ~ Chisq(n_tot - g)

#(delta_bar_i,j[k] - truedelta)/sqrt(S_pooled(1/n_i + 1/n_j)[k,k]) ~ t(n_tot - g)

#So:
alphaB <- 0.1/6
n_1 + n_2 + n_3  - 3

q.stud <- qt(1-alphaB/2, n_1 + n_2 + n_3  - 3)

S1 <- cov(stoneflakes[which(clust==1),])
S2 <- cov(stoneflakes[which(clust==2),])
S3 <- cov(stoneflakes[which(clust==3),])
S_p <- ((n_1 - 1)*S1 + (n_2 - 1)*S2 + (n_3 - 1)*S3)/(n_1 + n_2 + n_3  - 3)
band <- sqrt(c(diag(S_p)*(1/n_1 + 1/n_2),diag(S_p)*(1/n_1 + 1/n_3), diag(S_p)*(1/n_2 + 1/n_3)))*q.stud
smean_1<- sapply(stoneflakes[which(clust==1),], mean)
smean_2<- sapply(stoneflakes[which(clust==2),], mean)
smean_3<- sapply(stoneflakes[which(clust==3),], mean)
centers <- c(smean_1[1] - smean_2[1], smean_1[2] - smean_2[2], smean_1[1] - smean_3[1], smean_1[2] - smean_3[2], smean_2[1] - smean_3[1], smean_2[2] - smean_3[2])
CI <- cbind(inf = centers - band,
            centers = centers,
            sup = centers + band)
rownames(CI) <- c("L1 - L2", "W1 - W2", "L1 - L3", "W1 - W3", "L2 - L3", "W2 - W3")
CI
#Group 1 vs group 2
#There is difference in just width, width of groups 2 is higher

#Group 1 vs group 3
#Very different, group 1 is on average longer but group 3 is on average wider

#Group 2 vs group 3
#Very different, group 2 stoneflakes are largerer than those of group 3
