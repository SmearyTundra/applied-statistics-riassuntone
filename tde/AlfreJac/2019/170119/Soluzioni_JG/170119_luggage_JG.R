luggage <- read.table("luggage.txt", header = T)
n <- dim(luggage)[1]
p <- dim(luggage)[2]


#we can already assume gaussianity

#X ~ N4(means, SIGMA)

#Test:
#H0: mean.M.out + mean.M.ret = mean.F.out + mean.F.ret vs H1
#using a t-test
S <- cov(luggage)
means <- sapply(luggage, mean)

vec1 <- c(1,-1,1,-1)/2
T_stat1 <- abs(means%*%vec1)/sqrt(vec1%*%S%*%vec1/n)
p_val1 <- 1 - pt(T_stat1, n-1)
p_val1

#so we have strong evidence to say that there is a difference in luggage weight between men and women

#Test:
#H0: mean.M.out + mean.F.out = mean.M.ret + mean.F.ret vs H1
#using a t-test
vec2 <- c(1,1,-1,-1)/2
T_stat2 <- abs(means%*%vec2)/sqrt(vec2%*%S%*%vec2/n)
p_val2 <- 1 - pt(T_stat2, n-1)
p_val2

#so we have strong evidence to say that there is a difference in luggage weight between returning and outbound flights

alpha <- 0.1
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha, p, n-p)

CI.1 <- cbind(inf = vec1%*%means - sqrt(cfr.fisher)*sqrt(vec1%*%S%*%vec1/n),
              center = vec1%*%means,
              sup = vec1%*%means + sqrt(cfr.fisher)*sqrt(vec1%*%S%*%vec1/n))
CI.2 <- cbind(inf = vec2%*%means - sqrt(cfr.fisher)*sqrt(vec2%*%S%*%vec2/n),
              center = vec2%*%means,
              sup = vec2%*%means + sqrt(cfr.fisher)*sqrt(vec2%*%S%*%vec2/n))
CI <- rbind(m_f = CI.1,
            o_r = CI.2)
colnames(CI) <- c("inf", "center", "sup")
rownames(CI) <- c("m_f", "o_r")
CI
#indeed neither of the two intervals contain zero
#moreover we can say that, on average, men check in lighter luggages than women, and returning trips luggages are, again on average, heavier

#it is a prediction interval, in order to do that we need to consider not only the uncertainty for the mean of "female.returning", but also the variability around said mean
#since it is a multivariate gaussian, each single variable is still gaussian, so we want to predict
# x ~ N(mu, sigma2)
# x - xbar ~ N(0, (1+1/n)sigma2)
# (x-xbar)/sqrt((1+1/n)sigma2) ~ N(0,1)
# (n-1)S2/sigma2 ~ Chi(n-1)

# so (x-xbar)/sqrt((1+1/n)S2) ~ t(n-1)

smean.pred <- means[4]
s <- S[4,4]
band <- qt(1-alpha/2, n-1)*sqrt((1+1/n)*s)
CI.pred <- cbind(inf = smean.pred - band,
                 center = smean.pred,
                 sup = smean.pred + band)
CI.pred
#even though the point estimate for the expected value for the weight of the luggage of said person is smaller than the limit,
#we cannot exclude, with probability 90%, that said person will not be charged extra